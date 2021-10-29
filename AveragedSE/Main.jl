push!(LOAD_PATH, "/home/gridsan/wyao/Research/TraceFormula/Module")
#push!(LOAD_PATH, "/Users/jayyao/Documents/Research/TraceFormula/Module")
using GridapEM
using Gridap
using DelimitedFiles
using KrylovKit
using LinearAlgebra

"""This part is used to define all parameters used"""
# Geometry parameters of the mesh
# Rectangular with center circle design domain
λ = 1.0           # Wavelength
L = 2.5           # Width of the rectangular domain
H = 2.5           # Height of the rectangular domain
rd = 0.7          # Radius of the design domain circle
rt = rd + 0.2     # Radius of the target circle
dpml = 0.5        # Thickness of PML

# Characteristic length (controls the resolution, smaller the finer)
resol = 40.0      # Number of points per wavelength
l1 = λ / resol    # Normal region
l2 = l1 / 2.0     # Design region
l3 = 2 * l1       # PML

# Physical parameters 
k = 2 * π / λ       # Wave number 
# Bloch wavevector
kb = VectorValue(0.0, 0.0)     
ω = k               # c=1
ϵ1 = 1.0            # Relative electric permittivity for material 1 (y > 0)
ϵ2 = 1.0            # Relative electric permittivity for material 2 (y < 0)
ϵ3 = 0.0            # Relative electric permittivity for potential waveguide
ϵd = 12.0           # Relative electric permittivity for design material
μ = 1.0             # Relative magnetic permeability for all materials
wg_center = [0, 0]  # Waveguide center if exist
wg_size = [0, 0]    # Waveguide size if exist
#LHp = [Inf, h1 + h2]
#LHn = [Inf, h3]

# PML parameters
R = 1e-10           # Tolerence for PML reflection 
σ1 = -3 / 4 * log(R) / dpml / √ϵ1
σ2 = -3 / 4 * log(R) / dpml / √ϵ2
σs = [σ1, σ2]
############  Optimization parameters #############
flag_f = true       # Turn on filter
flag_t = true       # Turn on threshold

# Filter and threshold paramters
r = [0.02 * λ, 0.02 * λ]  # Filter radius
β = 80.0                  # β∈[1,∞], threshold sharpness
η = 0.5                   # η∈[0,1], threshold center

α = 1.0 / (2 * 1000.0)    # Equivalent loss α = 1/2Q

# Number of subspace
K = 20

# Amplify g for NLopt
Amp = 1

# Sum over kx
nkx = 30
nparts = nkx / 2

Bρ = true          # Matrix B depend on parameters?
ρv = 0.5

# Foundary constraint parameters
c = 0#resol^4
lw = r[1]
ls = r[1]
ηe = fηe(lw / r[1])
ηd = fηd(lw / r[1])


# Create mesh file
geo_param = CirRecGeometry(L, H, rd, rt, dpml, l1, l2, l3)
#geo_param = PeriodicGeometry(L, h1, h2, h3, ht, hs, dpml, l1, l2, l3)
meshfile_name = "geometry.msh"
MeshGenerator(geo_param, meshfile_name)

# Apply gridap finite element analysis to mesh file
gridap = GridapFE(meshfile_name, 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], [], flag_f)
#run(`gmsh geometry.msh`)

# Change default physics parameters
LHp=[L / 2, H / 2]  # Start of PML for x,y > 0
LHn=[L / 2, H / 2]  # Start of PML for x,y < 0

phys = PhysicalParameters(k, kb, ω, ϵ1, ϵ2, ϵ3, ϵd, μ, R, σs, dpml, LHp, LHn, wg_center, wg_size)

control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bρ, ρv, c, ηe, ηd)


Nri = 101
α = 1.0 / (2 * 100)
control = ControllingParameters(flag_f, flag_t, r, β, η, α, nparts, nkx, K, Amp, Bρ, ρv, c, ηe, ηd)
Powers = zeros(Nri)
Rs = zeros(Nri)
Ncv = zeros(Nri)
Ncv2 = zeros(Nri)
for ri = 1 : Nri
    rd = (ri - 1) * 0.9 / (Nri - 1) + 0.1
    Rs[ri] = rd
    rt = rd + 0.2
    geo_param = CirRecGeometry(L, H, rd, rt, dpml, l1, l2, l3)
    meshfile_name = "geometry.msh"
    MeshGenerator(geo_param, meshfile_name)
    gridap = GridapFE(meshfile_name, 1, 2, ["DirichletEdges", "DirichletNodes"], ["DesignNodes", "DesignEdges"], ["Target"], [], flag_f)
    
    N = num_free_dofs(gridap.FE_U)
    ρ0 = ones(gridap.np)
    ρf_vec = ρf_ρ0(ρ0; control, gridap)
    ρfh = FEFunction(gridap.FE_Pf, ρf_vec)
    ρth = (ρf -> Threshold(ρf; control)) ∘ ρfh
        
    A_mat = MatrixA(ρth; phys, control, gridap)
    B_mat = MatrixB(ρth; control, gridap)
    #@show sum(∫(ρth)gridap.dΩ_d) / sum(∫(1)gridap.dΩ_d)

    O_mat = MatrixOc(phys.k, phys.ϵ1; gridap)
    
    Neig = Int(ceil(ri / Nri * 30))
    G_ii, W_raw, info = eigsolve(x -> MatrixG(x; A_mat, B_mat, O_mat), rand(ComplexF64, N), Neig, :LM; krylovdim = 50)
    Ncv[ri] = num_contributing_values(G_ii, 0.99)
    Ncv2[ri] = num_contributing_values(G_ii, 0.9)
    @show Powers[ri] = sum(abs.(G_ii))
end