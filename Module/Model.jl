"""This file defines all models"""
struct PhysicalParameters
    k::Float64                   # Wavenumber of free space
    kb::VectorValue{2,Float64}   # Bloch wavevector
    ω::Float64                   # Angular frequency
    ϵ1::Float64                  # Electric permittivity for material 1 y > 0
    ϵ2::Float64                  # Electric permittivity for material 2 y < 0
    ϵ3::Float64                  # Electric permittivity for potential waveguide
    ϵd::Float64                  # Electric permittivity for design material
    μ::Float64                   # Magnetic permeability
    R::Float64                   # Reflection of PML
    σs::Vector{Float64}          # PML parameter
    dpml::Float64                # Thickness of PML
    LHp::Vector{Float64}         # Start of PML for x, y > 0
    LHn::Vector{Float64}         # Start of PML for x, y < 0
    wg_center::Vector{Float64}   # Center of waveguide
    wg_size::Vector{Float64}     # Size of waveguide
end

# PML coordinate streching functions
function s_PML(x; phys)
    σ = x[2]>0 ? phys.σs[1] : phys.σs[2]
    xf = [x[1], x[2]]
    u = @. ifelse(xf > 0 , xf - phys.LHp, - xf - phys.LHn)
    return @. ifelse(u > 0,  1 + (1im * σ / phys.k) * (u / phys.dpml)^2, $(1.0+0im))
end

function ds_PML(x; phys)
    σ = x[2]>0 ? phys.σs[1] : phys.σs[2]
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

function ξ0(x; phys)
    if phys.ϵ3 == 0
        return x[2] > 0 ? 1 / phys.ϵ1 : 1 / phys.ϵ2
    elseif abs(x[1] - phys.wg_center[1]) <= phys.wg_size[1] / 2 && abs(x[2] - phys.wg_center[2]) <= phys.wg_size[2] / 2
        return 1 / phys.ϵ3
    else
        return 1 / phys.ϵ1
    end
end

ξd(ρ, ϵmin, ϵmax, α)= 1 / ((ϵmin + (ϵmax - ϵmin) * ρ) * (1 + α * 1im)) - 1/ϵmin # in the design region

a_base(u, v; phys) = (x -> ξ0(x; phys)) * ((∇ .* (Λ(phys) * v) - 1im * phys.kb * v) ⊙ ((Λ(phys) .* ∇(u) + 1im * phys.kb * u))) - (phys.k^2 * phys.μ * (v * u))

a_design(u, v, ρth; phys, control) = ((ρ -> ξd(ρ, phys.ϵ1, phys.ϵd, control.α)) ∘ ρth) * ((∇(v) - 1im *phys.kb * v) ⊙ (∇(u) + 1im *phys.kb * u)) - phys.k^2 * 1im * control.α * phys.μ * (v * u)

a_gram(u, v; phys) = (phys.k^2 * phys.μ * (v * u))

function MatrixA(ρth; phys, control, gridap)
    A_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫(a_base(u, v; phys))gridap.dΩ + ∫(a_design(u, v, ρth; phys, control))gridap.dΩ_d
    end
    return lu(A_mat)
end

function MatrixB(ρbh; control, gridap)
    if control.Bρ
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            ∫(ρbh * (∇(u) ⋅ ∇(v)))gridap.dΩ_d
        end
    else
        B_mat = assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
            ∫((∇(u) ⋅ ∇(v)))gridap.dΓ_s
        end
    end
    return B_mat
end

function MatrixA0(phys, control, gridap)
    return assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫(a_gram(u, v; phys))gridap.dΩ + ∫(phys.k^2 * 1im * control.α * phys.μ * (v * u))gridap.dΩ_d
    end
end


Sp(u, v, ω, ϵ) = 1im / (4 * ω * ϵ) * (u * ∇(v) - v * ∇(u))

# Objective matrix for total emitted power from a circle
function MatrixOc(ω, ϵ; gridap)
    # Assemble the matrix
    return assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫( Sp(u, v, ω, ϵ) ⋅ (x -> (x / norm(x))) )gridap.dΓ_t[1]
    end
end

# Objective matrix for emitted power from horizontal line
function MatrixOl(ω, ϵ; gridap)
    # Assemble the matrix
    return assemble_matrix(gridap.FE_U, gridap.FE_V) do u, v
        ∫( Sp(u, v, ω, ϵ) ⋅ VectorValue(0, 1) )gridap.dΓ_t[1]
    end
end

# Objective vector for a waveguide mode overlap
function VectorO(Ey_eig, Mode_norm; gridap)
    l_temp(v) = ∫(v * Ey_eig / Mode_norm)gridap.dΓ_t[1]
    o_vec = assemble_vector(l_temp, gridap.FE_V)
    return o_vec
end
