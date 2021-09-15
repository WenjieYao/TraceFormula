"""
This file contains filter, threshold, volume constraint and foundary contraint functions. 
"""
struct ControllingParameters
    flag_f::Bool             # Enable filter
    flag_t::Bool             # Enable threshold
    r::Vector{Float64}       # r = [rx, ry] filter radius
    β::Float64               # Threshold steepness
    η::Float64               # Threshold value
    α::Vector{Float64}       # Equivalent loss term
    nparts::Int64            # Number of parts for paralell computing
    nkx::Int64               # Number of k points for k integral
    K::Int64                 # Number of contributing eigenvalues
    Amp::Float64             # Tuning amplitude for optimization
    Dρ::Bool                 # Matrix D depend on design parameter or not
    ρv::Float64              # Volume constraint on design parameter ρ 
    c::Float64               # Foundary constraint parameter
    ηe::Float64              # Foundary constraint parameter
    ηd::Float64              # Foundary constraint parameter
end
# ρf_vec = Filter(ρ_vec)
#a_f(r,u,v) = ∇(v)⊙(TensorValue(r[1]^2,0,0,r[2]^2)⋅∇(u))+v⊙u
a_f(r, u, v) = ∇(v) ⊙ (TensorValue(r[1]^2, 0, 0, r[2]^2) ⋅ ∇(u))

function Filter(ρ_vec; control, gridap)
    if (control.flag_f)
        ρh = FEFunction(gridap.FE_P, ρ_vec)
        op = AffineFEOperator(gridap.FE_Pf, gridap.FE_Qf) do u, v
            ∫(a_f(control.r, u, v))gridap.dΩ_d + ∫(v * u)gridap.dΩ, ∫(v * ρh)gridap.dΩ#, ∫( 0*v )gridap.dΓ_d
          end
        ρfh = solve(op)
        return get_free_values(ρfh)
    else
        return ρ_vec
    end
end

# Threshold function
#Threshold(pf;control) = control.flag_t==false ? pf : ((tanh(control.β*control.η)+tanh(control.β*(pf-control.η)))/(tanh(control.β*control.η)+tanh(control.β*(1.0-control.η))))
function Threshold(ρfh; control)
    if control.flag_t
        return ((tanh(control.β * control.η) + tanh(control.β * (ρfh - control.η)))/
                (tanh(control.β * control.η) + tanh(control.β * (1.0 - control.η))))
    else
        return ρfh
    end
end


function VolumeConstraint(ρW::Vector, grad::Vector; control, gridap)
    ρ0 = ρW[1 : gridap.np]
    ρf_vec = ρ_filter(ρ0; control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    
    if length(grad) > 0
        grad[gridap.np + 1 : end] = zeros(length(ρW) - gridap.np)
        l_temp(v) = ∫(v * ((ρf -> dρtdρf(ρf; control)) ∘ ρfh))gridap.dΩ_d
        grad0 = assemble_vector(l_temp, gridap.FE_Pf)
        grad[1 : gridap.np] = Dgdρ(grad0; control, gridap)
    end
    
    ρh = (ρf -> Threshold(ρf; control)) ∘ ρfh
    sum(∫(ρh)gridap.dΩ_d) - sum(∫(control.ρv)gridap.dΩ_d)
end

function fηe(x)
    if (x >= 0) && (x <= 1)
        return  0.25 * x^2 + 0.5
    elseif (x > 1) && (x <= 2)
        return -0.25 * x^2 + 1.0
    else
        return 1.0
    end
end

function fηd(x)
    if (x >= 0) && (x <= 1)
        return -0.25 * x^2 + 0.5
    elseif (x > 1) && (x <= 2)
        return  0.25 * x^2 + 1.0 - x
    else
        return 0
    end
end

fg(g, c) = exp(- c * norm(g)^2)
gc_LW(ph, ηe) = ph > ηe ? 0.0 : (ph - ηe)^2
gc_LS(ph, ηd) = ph < ηd ? 0.0 : (ph - ηd)^2

function LWConstraint(ρW::Vector, grad::Vector; control, gridap)
    ρ0 = ρW[1 : gridap.np]
    ρf_vec = ρ_filter(ρ0; control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρh = (ρf -> Threshold(ρf; control)) ∘ ρfh
    if length(grad) > 0
        grad[gridap.np + 1 : end] = zeros(length(ρW) - gridap.np)
        l_temp(v) = ∫(v * ((ph -> gc_LW(ph, control.ηe)) ∘ ρfh) 
                    * ((g -> fg(g, control.c)) ∘ ∇(ρfh)) 
                    * (((ρf -> dρtdρf(ρf; control)) ∘ ρfh) 
                    + 2 * ρh / (ρfh - control.ηe)) 
                    - 2 * control.c * ρh * ((ph -> gc_LW(ph, control.ηe)) ∘ ρfh) 
                    * ((g -> fg(g, control.c)) ∘ ∇(ρfh)) * (∇(v) ⋅ ∇(ρfh)))gridap.dΩ_d
        grad0 = assemble_vector(l_temp, gridap.FE_Pf)
        grad[1 : gridap.np] = Dgdρ(grad0; control, gridap)
    end
    
    
    sum(∫(ρh * ((g -> fg(g, control.c)) ∘ ∇(ρfh)) * ((ph -> gc_LW(ph, control.ηe)) ∘ ρfh))gridap.dΩ_d)
end