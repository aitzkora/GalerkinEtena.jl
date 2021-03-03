struct Mesh1D
   points::Array{Float64}
   cells::Array{Array{Int64,1},1}
   function Mesh1D(a::Float64, b::Float64, N::Int64)
       new(collect(LinRange(a,b,N+1)), map(collect,collect(zip(1:N, 2:N+1))))
   end
   function Mesh1D(pts::Array{Float64,1}, c::Array{Array{Int64,1},1})
       new(pts, c)
   end
end


struct SimplexMesh{D} 
  dim::Int64
  points::Array{Float64,2}
  cells::Array{Int64,2}
  function SimplexMesh{D}(points::Array{Float64,2}, cells::Array{Int64,2}) where {D}
    @assert Val(D) isa Union{map(x->Val{x},1:3)...}
    @assert D == size(points, 2)
    @assert D+1 == size(cells, 2)
    new{D}(D, points, cells)
  end
end
