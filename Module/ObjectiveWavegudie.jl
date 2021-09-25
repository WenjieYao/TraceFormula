function VectorQ(Ey_eig, Mode_norm, gridap)
    l_temp(v) = ∫(v * Ey_eig / Mode_norm) * gridap.dΓ_t[1]
    q_vec = assemble_vector(l_temp, gridap.FE_V)
    return q_vec
end

function g_v(v_vec; B_mat)
    return real(v_vec' * B_mat * v_vec)
end

#v_vec = v_pf(pf)
function v_ρf(ρf_vec; q_vec, phys, control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
    A_mat = MatrixA(ρth; phys, control, gridap)
    v_vec = A_mat' \ q_vec
    v_vec
end

# Chain Rule : dg/dp = dg/dg*dg/dv*dv/dpf*dpf/dp
# dg/dv=dg/dg*dg/dv
function rrule(::typeof(g_v),v_vec; B_mat)
  function g_pullback(dgdg)
    NO_FIELDS, dgdg * (B_mat * v_vec)
  end
  g_v(v_vec; B_mat), g_pullback
end


# dg/dpf=dg/dv*dv/dpf
function rrule(::typeof(v_ρf), ρf_vec; q_vec, phys, control, gridap)
    v_vec = v_ρf(ρf_vec; q_vec, phys, control, gridap)
    function U_pullback(dgdv)
      NO_FIELDS, Dgvdρf(dgdv, v_vec, ρf_vec; phys, control, gridap)
    end
    v_vec, U_pullback
end
  
dDpWG(pfh,v2h,dp;control) = control.Dp*real(((pf->dptdpf(pf;control))∘pfh)*v2h*dp)


function Dgvdρf(dgdv, v_vec, ρf_vec; phys, control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
    A_mat = MatrixA(ρth; phys, control, gridap)
    w_vec = A_mat \ dgdv
    
    vdh = FEFunction(gridap.FE_U, conj(v_vec))
    vh = FEFunction(gridap.FE_U, v_vec)
    wh = FEFunction(gridap.FE_V, w_vec)
    l_temp(dρ) = ∫(real(DBdρf(ρfh, vdh, vh; control) - 2 * DAdρf(ρfh, vdh, wh; phys, control)) * dρ)gridap.dΩ_d
    return assemble_vector(l_temp, gridap.FE_Pf)
end

# Final objective function
function gv_ρ(ρ0::Vector; q_vec, B_mat, phys, control, gridap)
    ρf_vec = ρf_ρ0(ρ0; control, gridap)
    v_vec = v_ρf(ρf_vec; q_vec, phys, control, gridap)
    g_v(v_vec; B_mat)
end

function gv_ρ(ρ0::Vector, grad::Vector; q_vec, phys, control, gridap)
    ρf_vec = ρf_ρ0(ρ0; control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh

    B_mat = MatrixB(ρth; control, gridap)
    if length(grad) > 0
        dgvdρ, = Zygote.gradient(ρ -> gv_ρ(ρ; q_vec, B_mat, phys, control, gridap), ρ0)
        grad[:] = dgvdρ * control.Amp
    end
    g_value = gv_ρ(ρ0; q_vec, B_mat, phys, control, gridap)

    #@show g_value
    open("gvalue.txt", "a") do io
        write(io, "$g_value \n")
    end
    # tc = readdlm("tcount.txt", Int64)[1]
    # open("PV/pvalue$(tc).txt", "w") do iop
    #     for i=1:gridap.np
    #         x_temp = ρ0[i]
    #         write(iop, "$x_temp \n")
    #     end
    # end
    # tc +=1
    # open("tcount.txt", "w") do iop
    #     write(iop, "$tc \n")
    # end

    return g_value * control.Amp
end


function gvρ_optimize(ρ_init, q_vec, TOL = 1e-4, MAX_ITER = 500; phys, control, gridap)
    ##################### Optimize #################
    opt = Opt(:LD_MMA, gridap.np)
    lb = zeros(gridap.np)
    ub = ones(gridap.np)
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.ftol_rel = TOL
    opt.maxeval = MAX_ITER
    opt.max_objective = (ρ0, grad) -> gv_ρ(ρ0, grad; q_vec, phys, control, gridap)
    if (length(ρ_init)==0)
        ρ_initial = readdlm("ρ_opt_value.txt", Float64)
        ρ_initial = ρ_initial[:]
    else
        ρ_initial = ρ_init[:]
    end
    if control.ρv < 1
        inequality_constraint!(opt, (x, g) -> VolumeConstraint(x, g; control, gridap), 1e-2)
    end
    if control.c > 0
        equality_constraint!(opt, (x, g) -> LWConstraint(x, g; control, gridap), 1e-8)
    end

    (g_opt, ρ_opt, ret) = optimize(opt, ρ_initial)
    @show numevals = opt.numevals # the number of function evaluations
    
    return g_opt / control.Amp, ρ_opt
end