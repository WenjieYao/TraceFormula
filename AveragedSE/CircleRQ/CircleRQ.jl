using Gmsh
using GridapGmsh
using Gridap
using DelimitedFiles
using KrylovKit
using LinearAlgebra
using Gridap.Fields
import Gmsh:gmsh

include("Mesh_R.jl")
include("GridapFE.jl")
include("Model.jl")

"""This part is used to define all parameters used"""
# Geometry parameters of the mesh
# Rectangular with center circle design domain
max_r = 2.0
start_r = 0.1
resol = 200.0      # Number of points per wavelength

Nri = 191
Q_list = [10, 50, 100, 1000]
Powers = zeros(Nri, length(Q_list) + 1)
Ncv1 = zeros(Int64, Nri, length(Q_list))
Ncv2 = zeros(Int64, Nri, length(Q_list))



λ = 1.0           # Wavelength
L = (max_r + 0.5) * 2           # Width of the rectangular domain
H = (max_r + 0.5) * 2           # Height of the rectangular domain
# rd = 0.7          # Radius of the design domain circle
rt = max_r + 0.2     # Radius of the target circle
dpml = 0.5        # Thickness of PML
# Characteristic length (controls the resolution, smaller the finer)
l1 = λ / resol    # Normal region
l2 = l1 / 1.0     # Design region
l3 = 2 * l1       # PML
# Physical parameters 
k = 2 * π / λ       # Wave number 
# Bloch wavevector
ω = k               # c=1
ϵ1 = 1.0            # Relative electric permittivity for material 1 (y > 0)
ϵd = 12.0           # Relative electric permittivity for design material
μ = 1.0             # Relative magnetic permeability for all materials
LHp=[L / 2, H / 2]  # Start of PML for x,y > 0
LHn=[L / 2, H / 2]  # Start of PML for x,y < 0

# PML parameters
R = 1e-10           # Tolerence for PML reflection 
σ = - 3 / 4 * log(R) / dpml / √ϵ1
############  Optimization parameters #############

# α = 1.0 / (2 * 1000.0)    # Equivalent loss α = 1/2Q
meshfile = "geometry.msh"
geo_param = RecGeometry(L, H, rt, dpml, l1, l2, l3)

MeshGenerator(geo_param, meshfile)
gridap = GridapFE(meshfile, 1, 2, ["DirichletEdges", "DirichletNodes"], ["Target"])
phys = PhysicalParameters(k, ω, ϵ1, ϵd, μ, R, σ, dpml, LHp, LHn)

N = num_free_dofs(gridap.FE_U)
O_mat = MatrixOc(phys.k, phys.ϵ1; gridap)

if isfile("ncv1.txt")
    run(`rm ncv1.txt`)
end

if isfile("ncv2.txt")
    run(`rm ncv2.txt`)
end

if isfile("power.txt")
    run(`rm power.txt`)
end

for ri = 1 : Nri
    rd = (ri - 1) * (max_r - start_r) / (Nri - 1) + start_r
    #rd = max_r
    Powers[ri, 1] = rd
    Neig = Int(ceil(ri / Nri * 25 * max_r))
    damp = 5 * rd / l1
    for qi = 1 : length(Q_list)
        α = 1.0 / (2 * Q_list[qi])
        A_mat = MatrixA(α, damp, rd; phys, gridap)
        B_mat = MatrixB(damp, rd; gridap)

        G_ii, W_raw, info = eigsolve(x -> MatrixG(x; A_mat, B_mat, O_mat), rand(ComplexF64, N), Neig, :LM; krylovdim = Neig + 5)
        Ncv1[ri, qi] = num_contributing_values(G_ii, 0.99)
        Ncv2[ri, qi] = num_contributing_values(G_ii, 0.9)
        Powers[ri, qi + 1] = sum(abs.(G_ii))
    end
    
    open("ncv1.txt", "a") do io
        writedlm(io, reshape(Ncv1[ri, :], (1, length(Q_list))))
    end

    open("ncv2.txt", "a") do io
        writedlm(io, reshape(Ncv1[ri, :], (1, length(Q_list))))
    end

    open("power.txt", "a") do io
        writedlm(io, reshape(Powers[ri, :], (1, length(Q_list)+1)))
    end
end


