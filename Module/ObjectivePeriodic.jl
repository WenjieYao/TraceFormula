function MatrixGk(x::Vector, ρth; O_mat, kx_s, phys, control, gridap)
    y=zeros(eltype(x),length(x))
    B_mat = MatrixB(ρth; control, gridap)
    for ki=1:length(kx_s)
        kx = kx_s[ki]
        kb = VectorValue(kx, 0.0)    
        physk = PhysicalParameters(phys.k, kb, phys.ω, phys.ϵ1, phys.ϵ2, phys.ϵ3, phys.ϵd, phys.μ, phys.R, phys.σs, phys.dpml, phys.LHp, phys.LHn, phys.wg_center, phys.wg_size)
        A_matk = MatrixA(ρth; phys=physk, control, gridap)
        y += A_matk' \ (O_mat * (A_matk \ (B_mat * x)))
    end
    return y/length(kx_s)
end

# Parallel sum k
function g_ρWk(ρW::Vector, kx; O_mat, phys, control, gridap)
    N = num_free_dofs(gridap.FE_U)
    @assert length(ρW) == (gridap.np + 2 * N * control.K)
    ρ0 = zeros(gridap.np)
    for i = 1 : gridap.np
        ρ0[i] = ρW[i]
    end
    W_mat = reinterpret(ComplexF64, reshape(ρW[gridap.np + 1 : end], (2 * N, control.K)))

    kb = VectorValue(kx, 0.0)    
    physk = PhysicalParameters(phys.k, kb, phys.ω, phys.ϵ1, phys.ϵ2, phys.ϵ3, phys.ϵd, phys.μ, phys.R, phys.σs, phys.dpml, phys.LHp, phys.LHn, phys.wg_center, phys.wg_size)
        
    ρf_vec = ρf_ρ0(ρ0; control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
    
    A_mat = MatrixA(ρth; phys=physk, control, gridap)
    B_mat = MatrixB(ρth; control, gridap)

    dgdρ, = Zygote.gradient(ρ -> g_ρ(ρ; O_mat, W_mat, phys=physk, control, gridap), ρ0)
        
    dgdW = reinterpret(Float64, DgdW(A_mat, W_mat, B_mat, O_mat))
    
    g_value = g_ρ(ρ0; O_mat, W_mat, phys=physk, control, gridap)

    g_grad = zeros(length(ρW) + 1)
    g_grad[1] = g_value
    g_grad[2 : gridap.np + 1] = dgdρ
    g_grad[gridap.np + 2 : end] = 2 * dgdW[:]
    
    return g_grad
end


function gρW_sumk(ρW::Vector, grad::Vector;ids, kx_s, dkx, O_mat, phys, control, gridap)
    gk = map_parts(ids) do myids
        mykxs = kx_s[myids]
        mygk = map(kx -> g_ρWk(ρW, kx; O_mat, phys, control, gridap), mykxs)
        return sum(mygk)
    end
    
    gvalue = 0.0
    igrad = 0
    for gki in gk
        igrad += 1
        if igrad == 1
            gvalue = sum(gki) * dkx / 2 / π * control.Amp
        elseif (length(grad)>0)
            grad[igrad-1] = sum(gki) * dkx / 2 / π * control.Amp
        end
    end
    #gvalue = sum(gk)*dkx/2/π*control.Amp
    open("gvalue.txt", "a") do io
        write(io, "$(gvalue/control.Amp) \n")
    end
    # tc = readdlm("tcount.txt", Int64)[1]
    # open("PV/pvalue$(tc).txt", "w") do iop
    #     for i=1:gridap.np
    #         x_temp = ρW[i]
    #         write(iop, "$x_temp \n")
    #     end
    # end
    # tc +=1
    # open("tcount.txt", "w") do iop
    #     write(iop, "$tc \n")
    # end
    return gvalue
end

function gρWk_optimize(ρ_init, L_local, TOL = 1e-4, MAX_ITER = 500; geo_param, phys, control)
    kx_ini = -π / geo_param.L
    dkx = 2 * π / geo_param.L / control.nkx
    kx_end = π / geo_param.L - dkx
    kx_s = range(kx_ini, kx_end;length = control.nkx)
    backend = SequentialBackend()
    parts = get_part_ids(backend, control.nparts)
    prange = PRange(parts, control.nkx)
    ids = map_parts(get_lid_to_gid, prange.partition)

    if (L_local != geo_param.L)
        MeshGenerator(geo_param, "geometry.msh")
    end
    gridap = GridapFE("geometry.msh", 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], ["Source"], control.flag_f)

    # Assemble matrices
    N = num_free_dofs(gridap.FE_U)
    O_mat = MatrixOl(phys.k, phys.ϵ1; gridap)
    
    ##################### Optimize #################
    opt = Opt(:LD_MMA, gridap.np + 2 * N * control.K)
    lb = zeros(gridap.np + 2 * N * control.K)
    lb[gridap.np + 1 : end] = - ones(2 * N * control.K) * Inf
    ub = ones(gridap.np + 2 * N * control.K)
    ub[gridap.np + 1 : end] = ones(2 * N * control.K) * Inf
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (ρW, grad) -> gρW_sumk(ρW, grad; ids, kx_s, dkx, O_mat, phys, control, gridap)
    if (length(ρ_init) == 0)
        ρW_initial = readdlm("ρW_opt_value.txt", Float64)
        ρW_initial = ρW_initial[:]
    else
        if (length(ρ_init) == 1)
            ρ_initial = ones(gridap.np) * ρ_init[1]
        else
            ρ_initial = ρ_init
        end
        ρW_initial = zeros(gridap.np + 2 * N * control.K)
        ρW_initial[1 : gridap.np] = ρ_initial[:]
        ρf_vec = ρf_ρ0(ρ_initial; control, gridap)
        ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
        ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh

        # A_mat = MatrixA(ρth; phys, control, gridap)
        # B_mat = MatrixB(ρth; control, gridap)
        
        G_ii, W_raw, info = eigsolve(x -> MatrixGk(x, ρth; O_mat, kx_s, phys, control, gridap), rand(ComplexF64, N), 2, :LM; krylovdim = 30)
        W_mat = rand(ComplexF64, N, control.K)
        W_mat = rand(ComplexF64, N, control.K)
        for ib = 1 : 2
            W_mat[:, ib] = W_raw[ib]
        end
        W_mat = reinterpret(Float64, W_mat)
        ρW_initial[gridap.np + 1 : end] = W_mat[:]
        @show abs(sum(G_ii))
    end
    if control.ρv < 1
        inequality_constraint!(opt, (x, g) -> VolumeConstraint(x, g; control, gridap), 1e-2)
    end
    if control.c > 0
        equality_constraint!(opt, (x, g) -> LWConstraint(x, g; control, gridap), 1e-8)
    end
    (g_opt, ρW_opt, ret) = optimize(opt, ρW_initial)
    @show numevals = opt.numevals # the number of function evaluations
    
    return g_opt / control.Amp, ρW_opt
end