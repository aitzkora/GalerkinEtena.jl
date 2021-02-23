module GalerkinEtena

using LinearAlgebra
using SpecialFunctions
using Revise
using Test
using SparseArrays



struct QuadratureFormula
    points::Array{Float64}
    weights::Array{Float64}
end

#1D and general funcs 
export JacobiGQ, JacobiGL, integrate, JacobiP,
       Legendre, computeElementaryMatrices, Mesh1D, genGrid, connect1D,
       computeMask, rk4, Advec1D, rhs1D, DGDiscretization, Maxwell1D

# 2D funcs
export rsToAb, WarpFactor, xyToRs, nodes2D, computeElementaryMatrices2D, 𝓔

include("mesh.jl")
include("integrate.jl")
include("assemble.jl")
include("assemble2D.jl")
include("maxwell1D.jl")
include("advec1D.jl")
include("rk4.jl")

end # module
