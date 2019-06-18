
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
vander(x, n)
computes the Vandermonde matrix V(x₁, x₂, ... , xₚ) defined by
```math
V_{ij} = [x_j^(i-1)]_{ij} \\, \\forall \\, i=1,\\cdots,p \\, j = 1,\\cdots,\\n
```
optional argument n enable us to choose the number of rows
"""
function vander(x::Array{Float64}, n = size(x,1))
    return x'.^(0:n-1)
end
let v_check = [ 1 1 1; 0.5 0.3 0.4; 0.25 0.09 0.16 ]
    @test vander([0.5, 0.3, 0.4]) ≈ v_check atol=1e-12
end

"""
GLT(N)
computes the Gauß-Lobatto-Чебышёв points defined by 
```math
T_k = -\\cos \\left( \\frac{k-1}{n-1}\\pi \\right)
```
"""
function GLT(n::Int)
    return -cos.((0:n-1)*π/(n-1))
end

"""
JacobiP
evaluate the Jacobi polynomial of type (α,β) > -1 (α+β ≢ -1) at points x for order N and returns P[1:length(xp))]
 Note   : They are normalized to be orthonormal.
taken from nodal-dg matlab code
"""
function JacobiP(x,α::Float64,β::Float64,N::Int)
xp = copy(x)
# Turn points into row if needed.
PL = zeros(N+1,size(xp,1))
# alias
γ = gamma
# Initial values P_0(x) and P_1(x)
γ0 = 2. .^(α+β+1.)/(α+β+1)*γ(α+1)*γ(β+1)/γ(α+β+1)
PL[1,:] .= 1.0/√γ0
if (N==0)
   return PL
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
