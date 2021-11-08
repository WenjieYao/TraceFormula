"""This file defines all models"""
struct PhysicalParameters
    k::Float64                   # Wavenumber of free space
    ω::Float64                   # Angular frequency
    ϵ1::Float64                  # Electric permittivity for material 1 y > 0
    ϵd::Float64                  # Electric permittivity for design material
    μ::Float64                   # Magnetic permeability
    R::Float64                   # Reflection of PML
    σ::Float64                   # PML parameter
    dpml::Float64                # Thickness of PML
    LHp::Vector{Float64}         # Start of PML for x, y > 0
    LHn::Vector{Float64}         # Start of PML for x, y < 0
end

# PML coordinate streching functions
function s_PML(x; phys)
    σ = phys.σ
    xf = [x[1], x[2]]
    u = @. ifelse(xf > 0 , xf - phys.LHp, - xf - phys.LHn)
    return @. ifelse(u > 0,  1 + (1im * σ / phys.k) * (u / phys.dpml)^2, $(1.0+0im))
end

function ds_PML(x; phys)
    σ = phys.σ
    xf = [x[1], x[2]]
    u = @. ifelse(xf > 0 , xf - phys.LHp, - xf - phys.LHn)
    ds = @. ifelse(u > 0, (2im * σ / phys.k) * (1 / phys.dpml)^2 * u, $(0.0+0im))
    return ds.*sign.(xf)
end

struct Λ <: Function
    phys
end

function (Λf::Λ)(x)
    s_x,s_y = s_PML(x; Λf.phys)
    return VectorValue(1/s_x, 1/s_y)
end

Fields.∇(Λf::Λ) = x -> TensorValue{2, 2, ComplexF64}(-(Λf(x)[1])^2 * ds_PML(x; Λf.phys)[1], 0, 0, -(Λf(x)[2])^2 * ds_PML(x; Λf.phys)[2])

function circle_damped(x, damp, rd)
    r = sqrt(x[1]^2 + x[2]^2) / rd - 1
    return 1 - 1 / (1 + exp(-damp * r))

end
function ξ0(x, α, damp, rd; phys)
    if circle_damped(x, damp, rd) > 0.5
        return 1 / (phys.ϵ1 + (phys.ϵd - phys.ϵ1) * circle_damped(x, damp, rd)) / (1 + α * 1im)
    else
        return 1 / (phys.ϵ1 + (phys.ϵd - phys.ϵ1) * circle_damped(x, damp, rd)) + 0im
    end
end

a(u, v, α, damp, rd; phys) = (x -> ξ0(x, α, damp, rd; phys)) * ((∇ .* (Λ(phys) * v)) ⊙ ((Λ(phys) .* ∇(u)))) - (phys.k^2 * phys.μ * (v * u))

function MatrixA(α, damp, rd; phys, gridap)
    A_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫(a(u, v, α, damp, rd; phys))gridap.dΩ 
    end
    return lu(A_mat)
end

function MatrixB(damp, rd; gridap)
    B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫((x->circle_damped(x, damp, rd)) * (∇(u) ⋅ ∇(v)))gridap.dΩ
    end
    return B_mat
end



Sp(u, v, ω, ϵ) = 1im / (4 * ω * ϵ) * (u * ∇(v) - v * ∇(u))

# Objective matrix for total emitted power from a circle
function MatrixOc(ω, ϵ; gridap)
    # Assemble the matrix
    return assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫( Sp(u, v, ω, ϵ) ⋅ (x -> (x / norm(x))) )gridap.dΓ_t
    end
end


function MatrixG(x::Vector; A_mat, B_mat, O_mat)
    A_mat' \ (O_mat * (A_mat \ (B_mat * x)))
end

function num_contributing_values(Gvec::Vector, cutoff = 0.99)
    nmv = length(Gvec)
    gsum = sum(abs.(Gvec))
    gtemp = 0
    for i = 1 : length(Gvec)
        gtemp += abs(Gvec[i])
        if (gtemp) > cutoff * gsum
            nmv = i
            break
        end
    end
    return nmv
end