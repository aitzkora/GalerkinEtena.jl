
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

let q=JacobiGQ(0.,0.,10)
    @test integrate(x->x*x, q) ≈ 2/3. atol=1.e-12
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
