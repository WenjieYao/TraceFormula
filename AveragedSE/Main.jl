push!(LOAD_PATH, "/home/jayyao/TraceFormula/Module")
using GridapEM
using Gridap
using DelimitedFiles
using KrylovKit
using LinearAlgebra

include("/home/jayyao/TraceFormula/Module/DefaultParameters.jl")

# Change default geometry parameters
rd = 0.5 #/ sqrt(2)
rt = rd + 0.2
resol = 40.0      # Number of points per wavelength
l1 = λ/resol      # Normal region
l2 = l1/2.0       # Design region
l3 = 2*l1         # PML

# Create mesh file
geo_param = CirRecGeometry(L, H, rd, rt, dpml, l1, l2, l3)
#geo_param = PeriodicGeometry(L, h1, h2, h3, ht, hs, dpml, l1, l2, l3)
meshfile_name = "geometry.msh"
MeshGenerator(geo_param, meshfile_name)

# Apply gridap finite element analysis to mesh file
gridap = GridapFE(meshfile_name, 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], [], flag_f)
#run(`gmsh geometry.msh`)

# Change default physics parameters
kb = VectorValue(2*π*0.,0)
LHp=[L / 2, H / 2]  # Start of PML for x,y > 0
LHn=[L / 2, H / 2]  # Start of PML for x,y < 0

phys = PhysicalParameters(k, kb, ω, ϵ1, ϵ2, ϵ3, ϵd, μ, R, σs, dpml, LHp, LHn, wg_center, wg_size)

# Change default control parameters
Bρ = true
ρv = 0.5
β = 80.0

control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bρ, ρv, c, ηe, ηd)

ρ_circ(x, r) = (x[1]^2 + x[2]^2) < r^2 ? 1 : 0
r_init = (0.5 - 0.5 / sqrt(2)) * 0.2 + 0.5 / sqrt(2)
lc_temp(v) = ∫(v * x->ρ_circ(x, r_init))gridap.dΩ
ρc_vec = assemble_vector(lc_temp, gridap.FE_P)
ρ_init = ρ_extract(ρc_vec; gridap)
#ρ_init[ρ_init .< 0.5] .= 0
ρ_init[ρ_init .> 0] .= 0.5
@show sum(ρ_init) / gridap.np, maximum(ρ_init)

r = [0.02 * λ, 0.02 * λ]  # Filter radius
Q_list = [50, 100, 500, 1000, 1000, 1000, 1000]
#α_list = [1000, 1000, 1000, 1000, 1000]
β_list = [5.0, 10.0, 20.0, 40.0, 80.0, 80.0, 80.0]

g_opt = 0
for bi = 1 : 7
    β = β_list[bi]
    α = 1.0 / (2 * Q_list[bi])
    phys = PhysicalParameters(k, kb, ω, ϵ1, ϵ2, ϵ3, ϵd, μ, R, σs, dpml, LHp, LHn, wg_center, wg_size)
    control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bρ, ρv, c, ηe, ηd)

    if bi == 1
        g_opt, ρW_opt = gρW_optimize([], 1e-6, 200; phys,control, gridap)
    else
        g_opt, ρW_opt = gρW_optimize([], 1e-6, 200; phys, control, gridap)
    end
    if isfile("ρW_opt.value.txt")
        run(`rm ρW_opt_value.txt`)
    end
    open("ρW_opt_value.txt", "w") do iop
        for i = 1 : length(ρW_opt)
            ρW_temp = ρW_opt[i]
            write(iop, "$ρW_temp \n")
        end
    end
    open("g_opt_value.txt", "a") do io
        write(io, "$g_opt \n")
    end
end
@show g_opt