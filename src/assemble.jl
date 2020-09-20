"""
computeElementaryMatrices(ξ::Array{Float64,1}, Np::Int)
computes the elementary matrices  𝓥, 𝓓ᵣ on the Gauß-Lobatto grid on [-1,1]
"""
function computeElementaryMatrices(ξ::Array{Float64,1},N::Int)
    𝓥, 𝓥ᵣ = Legendre(ξ, N)
    𝓓ᵣ = 𝓥ᵣ / 𝓥
    return 𝓥, 𝓓ᵣ
end

"""
genGrid(m::Mesh1D, ξ::Array{Float64})
generate a (#ξ,K) matrix corresponding to all degree of freedom points :
```math
G_{ij} = x^k_i \\in D^k
```
"""
function genGrid(m::Mesh1D, ξ::Array{Float64})
    np = size(ξ, 1)
    va = (x->x[1]).(m.cells)
    vb = (x->x[2]).(m.cells)
    grid = ones(np, 1) .* m.points[va]' + (ξ .+ 1) ./ 2 .* (m.points[vb] - m.points[va])'
    return grid
end

"""
computeMask(ξ::Array{Float64,1})
retrieves the index of the boundary nodes on the reference element
"""
function computeMask(ξ::Array{Float64,1})
    nodePrecision = 1e-12
    m1 = findall(abs.(ξ .+ 1) .< nodePrecision)
    m2 = findall(abs.(ξ .- 1) .< nodePrecision)
    return [m1; m2]
end


"""
e2e, e2f = connect1D(m::Mesh1D)

constructs the Element to Element (e2e) and Element to Face (e2f) matrices :
e2eᵢⱼ = k ⇔ element i is connected to k trough its j face if k ≠ i
e2fᵢⱼ = l ⇔ element i is connected to k trough its j face corresponding to l face in element k
"""

function connect1D(m::Mesh1D)
    K = length(m.cells)
    Nfaces = length(m.cells[1]) # assume nb_faces = cst here 2
    Ntot = Nfaces * K
    F2V = sparse(collect(1:Ntot), vcat(m.cells...), ones(Int8, Ntot) )
    F2F = F2V*F2V' - spdiagm(0=>ones(Int8,Ntot)) # if Fᵢⱼ = 1, global face i is connected to global face j
    faces = findnz(F2F)
    # convert face global indices to (element,face) ordering
    i2s=CartesianIndices((1:Nfaces,1:K))
    f1 = i2s[faces[1]]
    f2 = i2s[faces[2]]
    s2i = LinearIndices((1:K,1:Nfaces))
    ind = [s2i[z[2],z[1]] for z in f1]
    # rearrange into (K,Nfaces) shaped arrays
    E2E = [1:K;]*ones(Int64, 1, Nfaces)
    E2F = ones(Int64, K, 1) * (1:Nfaces)'
    E2E[ind] = [z[2] for z in f2]
    E2F[ind] = [z[1] for z in f2]
    return  E2E, E2F
end


"""
[x, vmapM, vmapP] =  DGDiscretization(m::Mesh1D, ξ::Array{Float64,1})

return a DG discretization along the mesh and the local discretization ξ
x is of size (#ξ, K )
vmapP and vmapP is a vector of size 2*K, if we reshape into a (2,K), matrix, then
[vmapM[1, i], vmapM[2, i]] are the indices of internal vertices of the i element
[vmapP[1, i], vmapP[2, i]] are the indices of external vertices of the i element

"""
function DGDiscretization(m::Mesh1D, ξ::Array{Float64})
    Np = length(ξ)
    K = length(m.cells)
    Nfaces = length(m.cells[1])
    nodeIds = reshape(1:K*Np, Np, K)

    fmask = computeMask(ξ)
    vmapM = nodeIds[fmask[:],1:K]

    vmapP = zeros(Int64, Nfaces, K)

    e2e, e2f = connect1D(m)

    x = genGrid(m, ξ )

    for k₁ = 1: K
        for f₁ = 1: Nfaces
            k₂ = e2e[k₁, f₁]
            f₂ = e2f[k₁, f₁]
            # f₂ is the corresponding face of f₁ in the neighbour k₂ 
            idM = vmapM[f₁, k₁]
            idP = vmapM[f₂, k₂]
            if (norm(x[idM]-x[idP]) < 1e-12)
                vmapP[f₁, k₁] = idP
            end
        end
    end

    # flatten the matrices
    vmapP = vmapP[:]
    vmapM = vmapM[:]
    return  x, vmapM, vmapP
end
