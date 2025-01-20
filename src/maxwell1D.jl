"""
structure for modelizing the 1D Maxwell's equation
```math
\\varepsilon(x)\\frac{\\partialE}{\\partial t}  = - \\frac{\\partialH}{\\partial x}  \\mbox{ and }
\\mu(x)\\frac{\\partialH}{\\partial t}  = - \\frac{\\partialE}{\\partial x} 
```
"""
struct Maxwell1D
    m::SimplexMesh{1}
    ξ::RefGrid{1}
    x::Matrix{Float64}
    vmapM::Vector{Int64}
    vmapP::Vector{Int64}
    vmapB::Vector{Int64}
    mapB::Vector{Int64}
    nx::Matrix{Float64}
    rx::Matrix{Float64}
    𝓓ᵣ::Matrix{Float64}
    fScale::Matrix{Float64}
    lift::Matrix{Float64}
    ε::Matrix{Float64}
    μ::Matrix{Float64}

    function Maxwell1D(m::SimplexMesh{1}, ξ::RefGrid{1})
        nx = normals(m)
        x, vmapM, vmapP = discretize(m, ξ)

        mapB  = findall( vmapM .== vmapP )
        vmapB = vmapM[ vmapM .== vmapP ]

        fmask = mask(ξ)

        𝓥, 𝓓ᵣ = elementaryMatrices(ξ)
        J = 𝓓ᵣ * x
        rx = 1. ./ J
        fScale = 1. ./ J[fmask, :]
        lift = 𝓥 * 𝓥' * 𝓔(fmask, ξ)
       
        K_2 = convert(Int64, m.K / 2)
        ε₁ = [ ones(1, K_2) 2*ones(1, K_2) ]
        μ₁ = ones(1, m.K)

        ε = ones(ξ.Np, 1) * ε₁
        μ = ones(ξ.Np, 1) * μ₁ 

        # create the object
        new(m, ξ, x, vmapM, vmapP, vmapB, mapB, nx, rx, 𝓓ᵣ, fScale, lift, ε, μ)
    end
end


"""
compute the right hand side of the maxwell equation u = [E, H]
"""
function rhs1D(pb::Maxwell1D, u::Matrix{Float64}, t::Float64)
    m_u_2 = convert(Int64, floor(size(u, 1) / 2))
    E = u[1:m_u_2,:]
    H = u[m_u_2+1:end,:]

    Zimp = sqrt.(pb.μ ./pb.ε)

    dE = zeros(2, pb.m.K)
    dE[:] = E[pb.vmapM] - E[pb.vmapP]

    dH = zeros(2, pb.m.K)
    dH[:] = H[pb.vmapM] - H[pb.vmapP]

    Zimpm = zeros(2, pb.m.K)
    Zimpm[:] = Zimp[pb.vmapM]

    Zimpp = zeros(2, pb.m.K)
    Zimpp[:] = Zimp[pb.vmapP]

    Yimpm = zeros(2, pb.m.K)
    Yimpm[:] = 1. ./ Zimpm[:]

    Yimpp = zeros(2, pb.m.K)
    Yimpp[:] = 1. ./ Zimpp[:]

    # Homogenenous Boundary Conditions Ez=0

    Ebc = -E[pb.vmapB] 
    dE[pb.mapB] = E[pb.vmapB] - Ebc
    Hbc = H[pb.vmapB]
    dH[pb.mapB] = H[pb.vmapB] - Hbc

    # evaluate upwind fluxes 
    fluxE = 1. ./ (Zimpm + Zimpp ) .* (pb.nx .* Zimpp .* dH - dE)
    fluxH = 1. ./ (Yimpm + Yimpp ) .* (pb.nx .* Yimpp .* dE - dH)

    # compute right hand sides of the PDE's 
    rhsE = (-pb.rx .* (pb.𝓓ᵣ*H) + pb.lift * (pb.fScale .* fluxE)) ./ pb.ε
    rhsH = (-pb.rx .* (pb.𝓓ᵣ*E) + pb.lift * (pb.fScale .* fluxH)) ./ pb.μ
    
    rhs = [ rhsE; rhsH ]
    return rhs
end

