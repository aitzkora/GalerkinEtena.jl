"""
computeElementaryMatrices2D(r::Array{Float64,1},s::Array{Float64,1}, N::Int)
computes the elementary matrices  𝓥, 𝓓ᵣ, 𝓓ₛ on the r,s coordinates corresponding to the integration grid
"""

function computeElementaryMatrices2D(r::Array{Float64,1}, s::Array{Float64,1},N::Int)

    Np = (N+1)*(N+2)÷2

    a,b = rsToAb(r,s)

    @assert length(r) == length(s)

    𝓥  = zeros(length(r), Np)
    𝓥ᵣ = zeros(length(r), Np)
    𝓥ₛ = zeros(length(s), Np)

    sk = 1
    for i=0:N
        for j=0:N-i 

           fa  = JacobiP(a, 0., 0., i)
           gb  = JacobiP(b, 2i+1., 0., j)
           # ∂ᵣJacobi(r,α,β,i) = √(i*(i+α+β+1))*Jacobi(r,α+1,β+1,i-1)
           dfa = i > 0 ? √(i*(i+1))*JacobiP(a, 1., 1., i-1) : 0.0
           dgb = j > 0 ? √(j*(j+2i+2.))*JacobiP(b, 2i+2., 1., j-1) : 0.0
           dr = dfa.* gb 
           ds = dfa .* (gb .* (1/2*(1 .+ a)))
           tmp = dgb .* ((1 .- b)/2).^i


           if (i > 0)
               factor = (1/2 * (1 .- b)).^(i - 1)
               dr = dr .* factor
               ds = ds .* factor
               tmp = tmp - 0.5 * i * gb .* factor
           end 

           ds += fa .* tmp

           # normalize
           n_fact = 2^(i+1/2)
           dr *= n_fact
           ds *= n_fact

           𝓥[:, sk] = √2fa.*gb.*(1 .-b ).^i;
           𝓥ᵣ[:, sk] = dr
           𝓥ₛ[:, sk] = ds

           sk += 1
        end
    end

    𝓓ᵣ = 𝓥ᵣ / 𝓥
    𝓓ₛ = 𝓥ₛ / 𝓥
    return 𝓥, 𝓓ᵣ, 𝓓ₛ
end


