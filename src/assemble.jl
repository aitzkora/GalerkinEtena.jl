"""
    elementaryMatrices(ξ::RefGrid{1})

computes the elementary matrices `𝓥, 𝓓ᵣ` on the Gauß-Lobatto grid on `I₁ = [-1,1]`
"""
function elementaryMatrices(ξ::RefGrid{1})
    𝓥, 𝓥ᵣ = Legendre(ξ.N, ξ.r[:,1])
    𝓓ᵣ = 𝓥ᵣ / 𝓥
    return 𝓥, 𝓓ᵣ
end

"""
    genGrid(m::SimplexMesh{1}, ξ::RefGrid{1})

generates a `(|ξ.r|,K)` matrix corresponding to all degree of freedom points :
```math
G_{ij} = x^k_i \\in D^k
```
"""

function genGrid(m::SimplexMesh{1}, ξ::RefGrid{1})
    va = m.cells[:,1]
    vb = m.cells[:,2]
    G = ones(ξ.Np, 1) .* m.points[va]' + (ξ.r .+ 1) ./ 2 .* (m.points[vb] - m.points[va])'
end

"""
    mask(ξ::RefGrid{1})

retrieves the index of the boundary nodes on the reference element
"""
function mask(ξ::RefGrid{1})
    nodePrecision = 1e-12
    # remember in 1D : [-1,1] is the ref element
    m1 = findall(abs.(ξ.r[:,1] .+ 1) .< nodePrecision) # "equal" to -1
    m2 = findall(abs.(ξ.r[:,1] .- 1) .< nodePrecision) # "equal" to 1
    return [m1; m2]
end

function vFace(m::SimplexMesh{1})
  return [1,2]
end

function vFace(m::SimplexMesh{2})
  return [1 2; 2 3; 3 1]
end

function vFace(m::SimplexMesh{3})
  return  [1 2 3; 1 2 4 ; 2 3 4 ; 1 3 4]
end

"""
    e2e, e2f = connect(m::SimplexMesh{D})

constructs the Element to Element (e2e) and Element to Face (e2f) sparse matrices :
- e2eᵢⱼ = k ⇔ element i is connected to k trough its j face if k ≠ i
- e2fᵢⱼ = l ⇔ element i is connected to k trough its j face corresponding to l face in element k
"""

function connect(m::SimplexMesh{D}) where D
  K = size(m.cells, 1)
  NFaces = D + 1 
  Ntot = NFaces * K
  Nv = maximum(m.cells)
  K = size(m.cells, 1)
  NFaces = size(m.cells, 2)
  Ntot = NFaces * K
  F2V = spzeros(Ntot, Nv)
  vf = vFace(m)
  sk = 1
  for k=1:K
    for face=1:NFaces
      F2V[sk, m.cells[k, vf[face,:]]] .= 1
      sk += 1
    end
  end

  F2F = F2V*F2V' - D * spdiagm(0=>ones(UInt8,Ntot)) # if F2Fᵢⱼ = 2, then i and j are the same edge (or face)
  faces = findnz(F2F.==D)
  # convert face global indices to (element,face) ordering
  i2s=CartesianIndices((1:NFaces,1:K))
  f1 = i2s[faces[1]]
  f2 = i2s[faces[2]]
  s2i = LinearIndices((1:K,1:NFaces))
  ind = [s2i[z[2],z[1]] for z in f1]
  # rearrange into (K,NFaces) shaped arrays
  E2E = [1:K;]*ones(Int64, 1, NFaces)
  E2F = ones(Int64, K, 1) * (1:NFaces)'
  E2E[ind] = [z[2] for z in f2]
  E2F[ind] = [z[1] for z in f2]
  return  E2E, E2F
end


function tiConnect(m::SimplexMesh{D}) where D
  Nfaces = D + 1
  K = size(m.cells, 1)
  fnodes =  generateFaces(m)
end 



"""
[x, vmapM, vmapP] =  discretize(m::SimplexMesh{1}, ξ::RefGrid{1})

return a DG discretization along the mesh and the local discretization ξ
x is of size (#ξ, K )
vmapP and vmapP is a vector of size 2*K, if we reshape into a (2,K), matrix, then
[vmapM[1, i], vmapM[2, i]] are the indices of internal vertices of the i element
[vmapP[1, i], vmapP[2, i]] are the indices of external vertices of the i element

"""
function discretize(m::SimplexMesh{1}, ξ::RefGrid{1})
    Np = ξ.Np
    K = size(m.cells, 1)
    @assert size(m.cells, 2) == ξ.NFaces
    nodeIds = reshape(1:K*Np, Np, K)

    fmask = mask(ξ)
    vmapM = nodeIds[fmask[:], 1:K]

    vmapP = zeros(Int64, ξ.NFaces, K)

    e2e, e2f = connect(m)

    x = genGrid(m, ξ )

    for k₁ = 1: K
        for f₁ = 1: ξ.NFaces
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
