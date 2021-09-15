module GridapEM

using Gmsh
using Gridap
using GridapGmsh
using SparseArrays
using ChainRulesCore
using Zygote
using LinearAlgebra
using KrylovKit
using PartitionedArrays
using NLopt
using DelimitedFiles

import ChainRulesCore: rrule
import Gmsh: gmsh

# Mesh functions
export CirRecGeometry
export PeriodicGeometry
export RecRecGeometry
export MeshGenerator

# Control functions
export ControllingParameters
export Threshold
export Filter
export VolumeConstraint
export LWConstraint
export fηe
export fηd

# Helper functions
export ρ_extend
export ρ_extract
export GaussianD
export GaussianY
export num_contributing_values

# Gridap FE function
export GridapParameters
export GridapFE

include("Mesh_Periodic.jl")
include("Mesh_CR.jl")
include("Mesh_RR.jl")
include("Control.jl")
include("Helper.jl")
include("GridapFE.jl")

end