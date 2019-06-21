"""
computeElementaryMatrices(Np::Int)
computes the elementary matrices 𝓜 , 𝓢 on the Gauß-Lobatto grid on [-1,1]
"""
function computeElementaryMatrices(ξ::Array{Float64},N::Int)
    𝓥, 𝓥ᵣ = Legendre(ξ, N)
    𝓓ᵣ=𝓥ᵣ / 𝓥
    𝓜  = inv(𝓥 * 𝓥')
    𝓢 = 𝓜  * 𝓓ᵣ
    return 𝓜 , 𝓢
end


function genGrid(ξ::Array{Float64}, m::Mesh1D)
    np = size(ξ, 1)
    va = m.cells[1]
    vb = m.cells[2]
    grid = ones(np, 1) * ξ[va] + (ξ .+ 1) ./ 2 .* (ξ[vb] - ξ[va])
end
