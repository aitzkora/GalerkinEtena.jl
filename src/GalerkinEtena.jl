module GalerkinEtena

using LinearAlgebra
using SpecialFunctions
using Revise
using Test
using SparseArrays

#1D and general funcs 
export JacobiGQ, JacobiGL, integrate, JacobiP, Legendre, 
       elementaryMatrices, SimplexMesh, RefGrid, Mesh1D, genGrid, connect,
       mask, rk4, advec1D, rhs1D, discretize, Maxwell1D, vFace,
       Advec, refGrid1D, RefGrid

# 2D funcs
export rsToAb, WarpFactor, xyToRs, nodes2D, computeElementaryMatrices2D, 𝓔

include("RefGrid.jl")
include("mesh.jl")
include("integrate.jl")
include("assemble.jl")
include("assemble2D.jl")
#include("maxwell1D.jl")
include("advec1D.jl")
include("rk4.jl")

end # module
