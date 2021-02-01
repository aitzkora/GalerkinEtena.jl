"""
JacobiGQ(α::Float64,β::Float64,N::Int)

return the N'th Gauß quadrature points
"""
function JacobiGQ(α::Float64,β::Float64,N::Int)
    if (N == 0)
        return QuadratureFormula(-(α-β)/(α+β+2.), 2.)
    end
    J = zeros(N+1)
    h₁ = 2*(0:N) .+ α .+ β
    J = diagm(0 => -1/2*(α^2-β^2)./(h₁.+2)./h₁)
    J += diagm(1 => 2. ./(h₁[1:N].+2.) .* .√((1:N).*((1:N) .+ α .+ β) .* ((1:N) .+ α) .* ((1:N) .+ β) ./ (h₁[1:N] .+ 1.) ./ (h₁[1:N] .+ 3.)))
    if (α + β < 10 * eps(1.) )
        J[1,1] = 0.0
    end 
    J = J + J';
    F = eigen!(J)
    points = F.values
    weights = (F.vectors[1,:]).^2*2^(α+β+1)/(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(β+α+1)
    return QuadratureFormula(points, weights)
end

"""
integrate(f,q)

computes the numerical approximation of the integral of f on the interval [-1,1] using q as a Quadrature formula

"""
function integrate(f, q::QuadratureFormula)
    return sum(f.(q.points).*q.weights)
end

"""
JacobiGL(α::Float64, β::Float64, N::Int)

Compute the N'th order Gauß Lobatto quadrature formula points
"""
function JacobiGL(α::Float64, β::Float64, N::Int)
    x = zeros(N + 1, 1)
    if (N == 1)
        return  [-1.,1.]
    end
    gq= JacobiGQ(α + 1, β + 1, N - 2)
    x = [-1; gq.points; 1]
end

"""
JacobiP(x::Array{Float64,1},α::Float64,β::Float64,N::Int)

evaluates the Jacobi polynomial of type (α,β) > -1 (α+β ≢ -1) at points x for order N 
Note : the Jacobi polynomial is normalize by a factor γₙ = √(2/(2n+1))
adapted from nodal-dg matlab code [https://github.com/tcew/nodal-dg]
"""
function JacobiP(x::Array{Float64,1},α::Float64,β::Float64,N::Int)
xp = copy(x)
# Turn points into row if needed.
PL = zeros(N+1,size(xp,1))
# alias
γ = gamma
# Initial values P_0(x) and P_1(x)
γ0 = 2. .^(α+β+1.)/(α+β+1)*γ(α+1)*γ(β+1)/γ(α+β+1)
PL[1,:] .= 1.0/√γ0
if (N==0)
    return PL[N+1,:]
end
γ1 = (α+1)*(β+1)/(α+β+3)*γ0;
PL[2,:] = ((α+β+2)*xp/2 .+ (α-β)/2)/√γ1;
if (N==1)
  return PL[N+1,:]
end

# Repeat value in recurrence.
a₋ = 2/(2+α+β)*√((α+1)*(β+1)/(α+β+3))

# Forward recurrence using the symmetry of the recurrence.
for i=1:N-1
  h1 = 2*i+α+β
  a₊ = 2/(h1+2.)*√((i+1)*(i+1+α+β)*(i+1+α)*(i+1+β)/(h1+1)/(h1+3))
  b₊ = - (α^2-β^2)/h1/(h1+2)
  PL[i+2,:] = 1/a₊*((xp.-b₊).*PL[i+1,:] -a₋*PL[i,:])
  a₋ =a₊
end
P = PL[N+1,:]
end


"""
Legendre(x::Array{Float64}, n::Int) computes the matrices
```math
P_{ij} = P^j(x_i)
```
and 
```math
P'_{ij} = \\frac{dP^j}{dx}(x_i)
```
where
```math
P^n(x) = \\frac{1}{2^nn!}\\frac{d^n}{dx^n}\\left((x^2-1)^n\\right)
```
Note 
"""
function Legendre(x::Array{Float64}, n::Int)
    γ = 2 ./(2(0:n).+1.)'
    P = zeros(size(x, 1), n + 1)
    P¹ = zeros(size(x, 1), n + 1)
    P[:, 1] = ones(size(x))
    P[:, 2] = x
    P¹[:, 1] = zeros(size(x))
    P¹[:, 2] = ones(size(x))
    for i=3:n+1
        P[:, i] = (2*i-3.)/(i-1) .* x .* P[:,i-1] - (i-2.)/(i-1) .* P[:, i - 2]
        P¹[:, i] = (2*i-3.)/(i-1) .* (P[:,i-1] + x .* P¹[:, i - 1]) - (i-2.)/(i-1) .* P¹[:, i - 2]
    end
    return P./ .√γ,P¹ ./ .√γ
end


"""
lagrange(α::Array{Float64})
computes the n² coefficients of the lagrange basis polynomial
taken from  Accuracy and Stability of Numerical Algorithms by Nicholas Higham p. 417
"""

function lagrange(α::Array{Float64})
    n = size(α,1)
    a = zeros(n+1)
    w = zeros(n,n)
    a[1] = -α[1]
    a[2] = 1
    for k=2:n
        a[k+1] = 1
        for j=k:-1:2
            a[j] = a[j-1] - α[k]*a[j]
        end
        a[1] = -α[k] * a[1]
    end
    for i = 1:n
        w[i,n] = 1
        s = 1
        for j=n-1:-1:1
            w[i,j]  = a[j+1] + α[i] * w[i, j+1]
            s = α[i]*s .+ w[i,j]
        end
        w[i, :] = w[i, :]/s
    end
    return w
end

# 2D functions

"""
rsToAb(r::Array{Float64,1}, s::Array{Float64,1})
changes (r,s) coordinates to (a,b) coordinates
"""
function rsToAb(r::Array{Float64, 1}, s::Array{Float64, 1})
    a = 2 * (1. .+ r) ./ (1 .- s ) .- 1.
    b = copy(s)
    a[s.==1.] .= -1.
    return a, b
end

function Vander2D(N::Int64, r::Array{Float64,1}, s::Array{Float64,1})
  a,b = rstoab(r,s)
  hcat([√2 * JacobiP(a,0.,0.,i) .* JacobiP(b, 0., 2*i+1., j) .* (1 .- b).^i 
        for i=0:N for j=0:N-i]...)
end
"""
WarpFactor(N::Int64, rout::Array{Float64,1})
compute the warping function as defined p. 177 in Warburton-Hesthaven
"""
function WarpFactor(N::Int64, rout::Array{Float64,1})

# Compute LGL and equidistant node distribution
LGLr = JacobiGL(0.,0.,N)
req  = LinRange(-1.,1.,N+1)

# Compute V based on req
𝓥 = hcat([JacobiP([req;],0.,0.,i) for i=0:N]...) # TODO : check if we can use legendre here

# Evaluate Lagrange polynomial at rout
Pmat = hcat([JacobiP(rout,0.,0.,i) for i=0:N]...)
Lmat = 𝓥'\Pmat'

# Compute warp factor
warp = Lmat'*(LGLr - req)

# Scale factor
zerof = (abs.(rout) .< (1.0-1.e-10))
sf = 1.0 .- (zerof.*rout).^2
warp = warp./sf + warp.*(zerof.-1.)

end

""" 
function xyToRs(x::Array{Float64,1}, y::Array{Float64,1})
convert (x,y) coords in equilateral triangle to (r,s) coordinates 
standard triangle I = [(-1,-1), (1,-1), (-1,1)]
"""
function xyToRs(x::Array{Float64,1}, y::Array{Float64,1})

λ1 = (sqrt(3.0)*y+1.0)/3.0
λ2 = (-3.0*x - sqrt(3.0)*y + 2.0)/6.0
λ3 = ( 3.0*x - sqrt(3.0)*y + 2.0)/6.0

return  -λ2 + λ3 - λ1, -λ2 - λ3 + λ1

end

"""
nodes2D(N::Int64)
computes the interpolation nodes (x,y) nodes in equilateral triangle for polynomial of order N             
"""
function nodes2D(N::Int64)

    α_opt = [0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999, 
             1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258]
    # Set optimized parameter, α depending on order N 
    if (N<16) 
        α = α_opt[N]
    else 
        α = 5/3.
    end

    # total number of nodes
    Np = (N+1)*(N+2)÷2

    # Create equidistributed nodes on equilateral triangle
    λ1 = zeros(Np) 
    λ2 = zeros(Np) 
    λ3 = zeros(Np)

    sk = 1
    for n=1:N+1
        for m=1:N+2-n
            λ1[sk] = (n-1)/N 
            λ3[sk] = (m-1)/N
            sk = sk+1;
        end
    end

    λ2 = 1.0 .-λ1 .-λ3;
    x = -λ2+λ3; y = (-λ2-λ3+2*λ1)/sqrt(3.0)

    # Compute blending function at each node for each edge
    blend1 = 4λ2.* λ3
    blend2 = 4λ1.* λ3
    blend3 = 4λ1.* λ2

    # Amount of warp for each node, for each edge
    warpf1 = WarpFactor(N, λ3 - λ2) 
    warpf2 = WarpFactor(N, λ1 - λ3) 
    warpf3 = WarpFactor(N, λ2 - λ1)

    # Combine blend & warp
    warp1 = blend1.*warpf1.*(1. .+ (α * λ1).^2)
    warp2 = blend2.*warpf2.*(1. .+ (α * λ2).^2)
    warp3 = blend3.*warpf3.*(1. .+ (α * λ3).^2)

    # Accumulate deformations associated with each edge
    x = x + 1*warp1 + cos(2*π/3)*warp2 + cos(4*π/3)*warp3
    y = y + 0*warp1 + sin(2*π/3)*warp2 + sin(4*π/3)*warp3
    
    return x,y
end
