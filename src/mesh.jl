struct SimplexMesh{D} 
  dim::Int64
  points::Array{Float64,2}
  cells::Array{Int64,2}
  K::Int64
  function SimplexMesh{D}(points::Array{Float64,2}, cells::Array{Int64,2}) where {D}
    @assert Val(D) isa Union{map(x->Val{x},1:3)...}
    @assert D == size(points, 2)
    @assert D+1 == size(cells, 2)
    K = size(cells, 1)
    new{D}(D, points, cells, K)
  end
end

function SimplexMesh{1}(v::Vector{Float64}, elem::Matrix{Int64})
  SimplexMesh{1}(v[:,1:1], elem)
end

function Mesh1D(a::Float64, b::Float64, N::Int64)
  SimplexMesh{1}(collect(LinRange(a,b,N+1)), vcat(map(x->[x[:]...]',zip(1:N, 2:N+1))...))
end

function normals(m::SimplexMesh{1})
  nx = [-ones(1, m.K); ones(1, m.K)];
end
