"""
computeElementaryMatrices(Np::Int)
computes the elementary matrices 𝓜 , 𝓢 on the Gauß-Lobatto grid on [-1,1]
"""
function computeElementaryMatrices(ξ::Array{Float64},N::Int)
    𝓥, 𝓥ᵣ = Legendre(ξ, N)
    𝓓ᵣ = 𝓥ᵣ / 𝓥
    𝓜  = inv(𝓥 * 𝓥')
    return 𝓜 , 𝓓ᵣ
end


function genGrid(ξ::Array{Float64}, m::Mesh1D)
    np = size(ξ, 1)
    va = (x->x[1]).(m.cells)
    vb = (x->x[2]).(m.cells)
    grid = ones(1, np) .* m.points[va] + (ξ' .+ 1) ./ 2 .* (m.points[vb] - m.points[va])
    return grid
end

function faces2Vertices(m::Mesh1D)
    K = length(m.cells)
    Nfaces = length(m.cells[1]) # assume nb_faces = cst here 2
    Ntot = Nfaces * K
    F2V = sparse(collect(1:Ntot), vcat(m.cells...), ones(Int64, Ntot) )
    F2F = F2V*F2V' - spdiagm(0=>ones(Int64,Ntot))
    faces = findnz(F2F)
    i2s=CartesianIndices((1:Nfaces,1:K))
    f1 = i2s[faces[1]]
    f2 = i2s[faces[2]]
    s2i = LinearIndices((1:K,1:Nfaces))
    ind = [s2i[z[2],z[1]] for z in f1]
    E2E = [1:K;]*ones(Int, 1, Nfaces)
    E2F = ones(Int, K, 1) * (1:Nfaces)'
    E2E[ind] = [z[2] for z in f2]
    E2F[ind] = [z[1] for z in f2]
    return  E2E, E2F
end


