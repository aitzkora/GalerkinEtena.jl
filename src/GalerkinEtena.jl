module GalerkinEtena

using LinearAlgebra
using SpecialFunctions
using Revise
using Test


struct QuadratureFormula
    points::Array{Float64}
    weights::Array{Float64}
end

export vander, GLT, JacobiGQ, JacobiGL, integrate, JacobiP, lagrange, legendre, diff_legendre

include("integrate.jl")
include("mesh.jl")

end # module
