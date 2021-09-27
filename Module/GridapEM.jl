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
using Gridap.Geometry
using Gridap.Fields
using Gridap.CellData: Interpolable

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
export Interpolated_Initial_Guess

# Gridap FE function
export GridapParameters
export GridapFE

# Model functions
export PhysicalParameters
export MatrixA
export MatrixB
export MatrixA0
export MatrixOc
export MatrixOl
export VectorO
export Sp

# Objective functions
export ρf_ρ0
export g_ρf
export g_ρ
export g_ρW
export MatrixG
export MatrixGk
export gρW_optimize

export gvρ_optimize
export gρW_sumk
export gρWk_optimize
export gv_ρ

include("Mesh_Periodic.jl")
include("Mesh_CR.jl")
include("Mesh_RR.jl")
include("Control.jl")
include("Helper.jl")
include("GridapFE.jl")
include("Model.jl")
include("Objective.jl")
include("ObjectivePeriodic.jl")
include("ObjectiveWavegudie.jl")
end