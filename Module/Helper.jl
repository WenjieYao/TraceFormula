"""Helper File""" 
# Convert piece-wise constant p (design region) to pvec (whole domain)
function ρ_extend(ρ0; gridap)
    ρ_vec = zeros(num_free_dofs(gridap.FE_P))
    ρi = 0
    @assert length(gridap.tags) == num_free_dofs(gridap.FE_P)
    for i = 1 : length(gridap.tags)
        if gridap.tags[i] == gridap.design_tag
            ρi += 1
            ρ_vec[i] = ρ0[ρi]
        end
    end
    ρ_vec
end

# Extract the design region part from a whole vector
function ρ_extract(ρ_vec; gridap)
    ρ0 = zeros(eltype(ρ_vec), gridap.np)
    ρi = 0
    @assert length(ρ_vec) == length(gridap.tags)
    for i = 1 : length(gridap.tags)
        if gridap.tags[i] == gridap.design_tag
            ρi += 1
            ρ0[ρi] = ρ_vec[i]
        end
    end
    @assert gridap.np == ρi
    ρ0
end

# Gaussian Distribution function with center x0
function GaussianD(x, x0::AbstractArray, δ::AbstractArray)
    n = length(x)
    @assert (n == length(x0)) && (n == length(δ))
    δn = 1.0
    x_δ = 0.0
    for i = 1 : n
        δn *= √(2π) * δ[i]
        x_δ += ((x[i] - x0[i]) / δ[i])^2
    end
    1.0 / δn * exp(- x_δ / 2.0)
end

# Gaussian Distribution function with center x0 in only Y direction
function GaussianY(x, x0::AbstractArray, δ::AbstractArray)
    n = length(x)
    @assert (n == length(x0)) && (n == length(δ))
    δn = √(2π) * δ[2]
    x_δ = ((x[2] - x0[2]) / δ[2])^2
    return abs(x[1] - x0[1]) <= (δ[1]) ? 1.0 / δn * exp(- x_δ / 2.0) : 0.0
end


# Evaluate the number of contributing values (sum/total > cutoff) for a vector
function num_contributing_values(Gvec::Vector, cutoff = 0.9)
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


function Interpolated_Initial_Guess(gridap)
    model_guess = GmshDiscreteModel("InitialGuess/geometry.msh")
    gridap_guess = GridapFE("InitialGuess/geometry.msh", 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], [], true)
    ρW_temp = readdlm("InitialGuess/ρW_opt_value.txt", Float64)
    ρW_temp = ρW_temp[:]
    ρ_init_guess = ρW_temp[1 : gridap_guess.np]
    ρh_guess = FEFunction(gridap_guess.FE_P, ρ_extend(ρ_init_guess; gridap=gridap_guess))
    ρfh_init = interpolate(Interpolable(ρh_guess),gridap.FE_Pf)
    ρh_init = interpolate(Interpolable(ρfh_init),gridap.FE_P)
    ρ_init = ρ_extract(get_free_dof_values(ρh_init); gridap)
    return ρ_init
end

