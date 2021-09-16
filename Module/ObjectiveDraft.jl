NO_FIELDS = ZeroTangent()
# Objective trace 
function g_U(U_mat; O_mat, WBW)
    g_temp = tr((U_mat' * O_mat * U_mat) / WBW)
    @assert abs(imag(g_temp) / real(g_temp)) < 1e-6
    real(g_temp)
end

#U_mat = U_ρf(ρf)
function U_ρf(ρf_vec; B_mat, W_mat, phys, control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
    A_mat = MatrixA(ρth; phys, control, gridap)
    U_mat = A_mat \ (B_mat * W_mat)
    U_mat
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

# Chain Rule : dg/dp = dg/dg*dg/du*du/dpf*dpf/dp
# dg/du=dg/dg*dg/du
function rrule(::typeof(g_U), U_mat; O_mat, WBW)
  function g_pullback(dgdg)
    NO_FIELDS, dgdg * DgdU(U_mat; O_mat, WBW)
  end
  g_U(U_mat;O_mat, WBW), g_pullback
end

function DgdU(U_mat; O_mat, WBW)
    (O_mat * U_mat) / WBW
end


# dg/dρf=dg/du*du/dρf
function rrule(::typeof(U_ρf), ρf_vec; B_mat, W_mat, phys, control, gridap)
    U_mat = U_ρf(ρf_vec; B_mat, W_mat, phys, control, gridap)
    function U_pullback(dgdU)
      NO_FIELDS, Dgdρf(dgdU, U_mat, W_mat, B_mat, ρf_vec; phys, control, gridap)
    end
    U_mat, U_pullback
end

Dρtdρf(ρf; control) = control.flag_t ? control.β * (1.0 - tanh(control.β * (ρf - control.η))^2) / (tanh(control.β * control.η) + tanh(control.β * (1.0 - control.η))) : 1.0

Dξdρf(ρf, ϵmin, ϵmax; control)= (ϵmin - ϵmax) / (ϵmin + (ϵmax - ϵmin) * Threshold(ρf; control))^2 / (1 + control.α * 1im) * Dρtdρf(ρf; control)

DAdρf(ρfh, u, v; phys, control, gridap) = ((ρf -> Dξdρf(ρf, phys.ϵ1, phys.ϵd; control)) ∘ ρfh) * ((∇(v) - 1im *phys.kb * v) ⊙ (∇(u) + 1im *phys.kb * u))

DBdρf(ρfh, u, v; control, gridap) = ((ρf -> Dρtdρf(ρf; control)) ∘ ρfh) * (∇(v) ⊙ ∇(u))

function Dgdρf(dgdU, U_mat, W_mat, B_mat, ρf_vec; phys, control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
    A_mat = MatrixA(ρth; phys, control, gridap)
    Z_mat = A_mat'\dgdU
    Wr_mat = W_mat / (W_mat' * B_mat * W_mat)
    Wl_mat = W_mat * dgdU' * U_mat
    
    dgdρf = zeros(num_free_dofs(gridap.FE_Pf))
    for k_i = 1 : size(dgdU)[2]
        uh = FEFunction(gridap.FE_U, U_mat[:, k_i])
        zh = FEFunction(gridap.FE_V, conj(Z_mat[:, k_i]))
        
        l_temp1(dρ) = ∫(- 2 * real(DAdρf(ρfh, uh, zh; phys, control, gridap)) * dρ)gridap.dΩ_d
        dgdρf += assemble_vector(l_temp1, gridap.FE_Pf)
        if control.Bρ
            wh = FEFunction(gridap.FE_U, W_mat[:, k_i])
            wrh = FEFunction(gridap.FE_U, Wr_mat[:, k_i])
            wlh = FEFunction(gridap.FE_V, conj(Wl_mat[:, k_i]))
            l_temp2(dρ) = ∫(real(2 * DBdρf(ρfh, wh, zh; control, gridap) - DBdρf(ρfh, wrh, wlh; control, gridap)) * dρ)gridap.dΩ_d
            dgdρf += assemble_vector(l_temp2, gridap.FE_Pf)
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
function g_ρ(ρ0::Vector; O_mat, B_mat, W_mat, phys, control, gridap)
    ρf_vec = ρf_ρ0(ρ0; control, gridap)
    U_mat = U_ρf(ρf_vec; B_mat, W_mat, phys, control, gridap)
    WBW = W_mat' * B_mat * W_mat
    g_U(U_mat; O_mat, WBW)
end

function g_ρ(ρ0::Vector, grad::Vector; O_mat, B_mat, W_mat, phys, control, gridap)
    if length(grad) > 0
        dgdρ, = Zygote.gradient(ρ -> g_ρ(ρ; O_mat, B_mat, W_mat, phys, control, gridap), ρ0)
        grad[:] = dgdρ
    end
    g_value = g_ρ(ρ0::Vector; O_mat, B_mat, W_mat, phys, control, gridap)
    return g_value
end