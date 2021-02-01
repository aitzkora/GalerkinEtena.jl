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

export vander, JacobiGQ, JacobiGL, integrate, JacobiP, lagrange,
       Legendre, computeElementaryMatrices, Mesh1D, genGrid, connect1D,
       computeMask, rk4, Advec1D, rhs1D, DGDiscretization, Maxwell1D, rsToAb, WarpFactor, xyToRs

include("mesh.jl")
include("integrate.jl")
include("assemble.jl")
include("maxwell1D.jl")
include("advec1D.jl")
include("rk4.jl")

end # module
