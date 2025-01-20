"""
type for discretisation in r,s,t coordinates
"""
struct RefGrid{D}
  dim::Int64
  N::Int64
  Np::Int64 
  Nfp::Int64
  NFaces::Int64
  r::Matrix{Float64}
  function RefGrid{D}(N::Int64, Np::Int64, Nfp::Int64, NFaces::Int64, r::Matrix{Float64}) where {D}
    @assert Val(D) isa Union{map(x->Val{x},1:3)...}
    new{D}(D, N, Np, Nfp, NFaces, r)
  end
end

"""
    RefGrid1D(a::Float64, b::Float64, N::Int64)

initialize a 1D grid reference on [a,b] with a N order polynomial
"""
function refGrid1D(a::Float64, b::Float64, N::Int64)
  v = JacobiGL(0., 0., N)
  RefGrid{1}(N, N+1, 1, 2, v[:,1:1])
end
