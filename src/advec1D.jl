"""
structure for modelizing the 1D advection equation

```math
\\frac{\\partial u}{\\partial t} + a \\frac{\\partial u}{\\partial x} = 0
```
with
```math
u(x,0) = \\sin (x) \\\\
u(0,t) = -\\sin (a t)
```
"""
struct Advec{D}
    m::SimplexMesh{D}
    ξ::RefGrid{D}
    x::Array{Float64,2}
    vmapM::Array{Int64,1}
    vmapP::Array{Int64,1}
    nx::Array{Float64,2}
    rx::Array{Float64,2}
    𝓓ᵣ::Array{Float64,2}
    fScale::Array{Float64,2}
    lift::Array{Float64,2}


    function Advec{D}(m::SimplexMesh{D}, ξ::RefGrid{D}) where D
        nx = normals(m)
        x, vmapM, vmapP = Discretize(m, ξ)
        fmask = mask(ξ)
        # compute the metric and jacobian
        𝓥, 𝓓ᵣ = computeElementaryMatrices(ξ)
        J = 𝓓ᵣ * x
        rx = 1. ./ J
        fScale = 1. ./ J[fmask, :]
        lift = 𝓥 * 𝓥' * 𝓔(fmask, ξ)
        # create the object
        new(m, ξ, x, vmapM, vmapP, nx, rx, 𝓓ᵣ, fScale, lift)
    end
end


function advec1D(a, b, K, Np)
  m = Mesh1D(a, b, K)
  ξ = RefGrid
  ad = Advec{1}

    return ad
end
"""
    rhs1D(ad::Advec1D, u::Array{Float64,2}, t::Float64, a::Float64, α::Float64)

computes the right hand side of the advection problem
"""
function rhs(ad::Advec{1}, u::Array{Float64,2}, t::Float64, a::Float64, α::Float64)
    mapI = 1
    mapO = ad.ξ.K * 2
    vmapI = 1
    vmapO = ad.ξ.K * ad.ξ.Np

    # compute the numerical flux using jumps
    du = zeros(2, ad.ξ.K)
    du[:] = (u[ad.vmapM] - u[ad.vmapP]) .* (a * ad.nx[:] - (1 - α) * abs.(a .* ad.nx[:])) ./ 2.

    #impose boundary condition at x = 0
    uᵢₙ = -sin.(a .* t)
    du[mapI] = (u[vmapI] - uᵢₙ ) .* (a .* ad.nx[mapI] - ( 1 - α ) * abs.( a .* ad.nx[mapI])) ./ 2
    du[mapO] = 0

    rhs = -a * ad.rx .* (ad.𝓓ᵣ * u ) + ad.lift * (ad.fScale .* du)
    return rhs
end

