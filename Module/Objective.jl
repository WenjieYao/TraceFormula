NO_FIELDS = ZeroTangent()
# Objective trace 
function g_ρf(ρf_vec; O_mat, W_mat, phys, control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
    A_mat = MatrixA(ρth; phys, control, gridap)
    B_mat = MatrixB(ρth; control, gridap)
    U_mat = A_mat \ (B_mat * W_mat)
    g_temp = tr((U_mat' * O_mat * U_mat) / (W_mat' * B_mat * W_mat))
    @assert abs(imag(g_temp) / real(g_temp)) < 1e-6
    real(g_temp)
end

#ρf = ρf_ρ0(ρ0)
function ρf_ρ0(ρ0; control, gridap)
    ρ_vec = ρ_extend(ρ0; gridap)
    ρf_vec = Filter(ρ_vec; control, gridap)
    ρf_vec[ρf_vec .< 0] .= 0
    ρf_vec[ρf_vec .> 1] .= 1.0
    ρf_vec
end

function MatrixG(x::Vector; A_mat, B_mat, O_mat)
    A_mat' \ (O_mat * (A_mat \ (B_mat * x)))
end

# Chain Rule : 
# dg/dρf=dg/dg * dg/dρf
function rrule(::typeof(g_ρf), ρf_vec; O_mat, W_mat, phys, control, gridap)
    function U_pullback(dgdg)
      NO_FIELDS, dgdg * Dgdρf(ρf_vec, O_mat, W_mat; phys, control, gridap)
    end
    g_ρf(ρf_vec; O_mat, W_mat, phys, control, gridap), U_pullback
end

Dρtdρf(ρf; control) = control.flag_t ? control.β * (1.0 - tanh(control.β * (ρf - control.η))^2) / (tanh(control.β * control.η) + tanh(control.β * (1.0 - control.η))) : 1.0

Dξdρf(ρf, ϵmin, ϵmax; control)= (ϵmin - ϵmax) / (ϵmin + (ϵmax - ϵmin) * Threshold(ρf; control))^2 / (1 + control.α * 1im) * Dρtdρf(ρf; control)

DAdρf(ρfh, u, v; phys, control) = ((ρf -> Dξdρf(ρf, phys.ϵ1, phys.ϵd; control)) ∘ ρfh) * ((∇(v) - 1im *phys.kb * v) ⊙ (∇(u) + 1im *phys.kb * u))

DBdρf(ρfh, u, v; control) = ((ρf -> Dρtdρf(ρf; control)) ∘ ρfh) * (∇(v) ⊙ ∇(u))

function Dgdρf(ρf_vec, O_mat, W_mat; phys, control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
    A_mat = MatrixA(ρth; phys, control, gridap)
    B_mat = MatrixB(ρth; control, gridap)
    
    U_mat = A_mat \ (B_mat * W_mat)
    WBW = W_mat' * B_mat * W_mat
    dgdU = (O_mat * U_mat) / WBW
    Z_mat = A_mat' \ dgdU
    #Z_mat = A_mat' \ ((O_mat * U_mat) / WBW)
    Wr_mat = W_mat / WBW
    Wl_mat = W_mat * (dgdU' * U_mat)
    #Wl_mat = W_mat / WBW * (U_mat' * O_mat * U_mat)
    
    dgdρf = zeros(num_free_dofs(gridap.FE_Pf))
    for k_i = 1 : control.K
        uh = FEFunction(gridap.FE_U, U_mat[:, k_i])
        zh = FEFunction(gridap.FE_V, conj(Z_mat[:, k_i]))
        if control.Bρ
            wh = FEFunction(gridap.FE_U, W_mat[:, k_i])
            wrh = FEFunction(gridap.FE_U, Wr_mat[:, k_i])
            wlh = FEFunction(gridap.FE_V, conj(Wl_mat[:, k_i]))
            l_temp2(dρ) = ∫(real(2 * DBdρf(ρfh, wh, zh; control) - 2 * DAdρf(ρfh, uh, zh; phys, control) - DBdρf(ρfh, wrh, wlh; control)) * dρ)gridap.dΩ_d
            dgdρf += assemble_vector(l_temp2, gridap.FE_Pf)
        else
            l_temp1(dρ) = ∫(- 2 * real(DAdρf(ρfh, uh, zh; phys, control)) * dρ)gridap.dΩ_d
            dgdρf += assemble_vector(l_temp1, gridap.FE_Pf)
        end
    end
    return dgdρf
end
        
# dg/dρ=dg/dρf*dρf/dρ
function rrule(::typeof(ρf_ρ0), ρ0; control, gridap)
  function ρf_pullback(dgdρf)
    NO_FIELDS, Dgdρ(dgdρf; control, gridap)
  end
  ρf_ρ0(ρ0; control, gridap), ρf_pullback
end

function Dgdρ(dgdρf; control, gridap)
    if control.flag_f
        Af = assemble_matrix(gridap.FE_Pf, gridap.FE_Qf) do u, v
            ∫(a_f(control.r, u, v))gridap.dΩ_d + ∫(v * u)gridap.dΩ
        end
        λvec = Af' \ dgdρf
        λh = FEFunction(gridap.FE_Pf, λvec)
        l_temp(dρ) = ∫(λh * dρ)gridap.dΩ
        return ρ_extract(assemble_vector(l_temp, gridap.FE_P); gridap)
    else
        return ρ_extract(dgdρf; gridap)
    end
end

# Final objective function
function g_ρ(ρ0::Vector; O_mat, W_mat, phys, control, gridap)
    ρf_vec = ρf_ρ0(ρ0; control, gridap)
    g_ρf(ρf_vec; O_mat, W_mat, phys, control, gridap)
end

function g_ρ(ρ0::Vector, grad::Vector; O_mat, W_mat, phys, control, gridap)
    if length(grad) > 0
        dgdρ, = Zygote.gradient(ρ -> g_ρ(ρ; O_mat, W_mat, phys, control, gridap), ρ0)
        grad[:] = dgdρ
    end
    g_value = g_ρ(ρ0::Vector; O_mat, W_mat, phys, control, gridap)
    return g_value
end

# Adding W dependence
function DgdW(A, W, B, O)
    WBW = W' * B * W
    U = A \ (B * W)
    B' * (A' \ (O * (U / WBW))) - (B * W / WBW) * (U' * (O * (U / WBW)))
end

function g_ρW(ρW::Vector, grad::Vector; O_mat, phys, control, gridap)
    N = num_free_dofs(gridap.FE_U)
    @assert length(ρW) == (gridap.np + 2 * N * control.K)
    ρ0 = zeros(gridap.np)
    for i = 1 : gridap.np
        ρ0[i] = ρW[i]
    end
    W_mat = reinterpret(ComplexF64, reshape(ρW[gridap.np + 1 : end], (2 * N, control.K)))

    if length(grad) > 0
        ρf_vec = ρf_ρ0(ρ0; control, gridap)
        ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
        ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
        
        A_mat = MatrixA(ρth; phys, control, gridap)
        B_mat = MatrixB(ρth; control, gridap)
        
        dgdρ, = Zygote.gradient(ρ -> g_ρ(ρ; O_mat, W_mat, phys, control, gridap), ρ0)
        grad[1 : gridap.np] = dgdρ * control.Amp

        dgdW = reinterpret(Float64, DgdW(A_mat, W_mat, B_mat, O_mat))
        grad[gridap.np + 1 : end] = 2 * control.Amp * dgdW[:]
    end
    g_value = g_ρ(ρ0; O_mat, W_mat, phys, control, gridap)
    open("gvalue.txt", "a") do io
        write(io, "$g_value \n")
    end
    return g_value * control.Amp
end


function gρW_optimize(ρ_init, TOL = 1e-4, MAX_ITER = 500, OptAlg = :LD_MMA, IsQl = false; phys, control, gridap)
    # Assemble matrices
    N = num_free_dofs(gridap.FE_U)
    if IsQl
        O_mat = MatrixOl(phys.k, phys.ϵ1; gridap)
    else
        O_mat = MatrixOc(phys.k, phys.ϵ1; gridap)
    end
    
    ##################### Optimize #################
    opt = Opt(OptAlg, gridap.np + 2 * N * control.K)
    lb = zeros(gridap.np + 2 * N * control.K)
    lb[gridap.np + 1 : end] = - ones(2 * N * control.K) * Inf
    ub = ones(gridap.np + 2 * N * control.K)
    ub[gridap.np + 1 : end] = ones(2 * N * control.K) * Inf
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (ρW, grad) -> g_ρW(ρW, grad; O_mat, phys, control, gridap)
    if (length(ρ_init) == 0)
        ρW_initial = readdlm("pW_opt_value.txt", Float64)
        ρW_initial = ρW_initial[:]
    else
        ρ_initial = ρ_init
        ρW_initial = zeros(gridap.np + 2 * N * control.K)
        ρW_initial[1 : gridap.np] = ρ_initial[:]
        ρf_vec = ρf_ρ0(ρ_initial; control, gridap)
        ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
        ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
        
        A_mat = MatrixA(ρth; phys, control, gridap)
        B_mat = MatrixB(ρth; control, gridap)
        
        G_ii, W_raw, info = eigsolve(x -> MatrixG(x; A_mat, B_mat, O_mat), rand(ComplexF64, N), min(5, control.K), :LM; krylovdim = 30)
        W_mat = rand(ComplexF64, N, control.K)
        for ib = 1 : min(5, control.K)
            W_mat[:, ib] = W_raw[ib]
        end
        W_mat = reinterpret(Float64, W_mat)
        ρW_initial[gridap.np + 1 : end] = W_mat[:]
    end
    if control.ρv < 1
        inequality_constraint!(opt, (x, g) -> VolumeConstraint(x, g; control, gridap), 1e-2)
    end
    if control.c > 0
        equality_constraint!(opt, (x, g) -> LWConstraint(x, g; control, gridap), 1e-8)
    end
    (g_opt, ρW_opt, ret) = optimize(opt, ρW_initial)
    
    return g_opt / control.Amp, ρW_opt
end