using GalerkinEtena
using LinearAlgebra
Np = 6
ξ=JacobiGL(0.,0.,Np-1)
q10=JacobiGQ(0.,0.,10)
𝓥 = hcat([JacobiP(ξ,0.,0.,i) for i=0:Np-1]...)
𝓜  = inv(𝓥*𝓥')
x=q10.points
ϕₖ=lagrange(ξ)
ϕ=(i,x)->(x'.^(0:Np-1))'*ϕₖ[i,:]
res=hcat([[sum(ϕ(i,x).*ϕ(j,x).*q10.weights) for i=1:Np] for j=1:Np]...)
println("|𝓜 -res| = ", norm(𝓜 -res))
r = [0.65]
Iᵣ = vcat([ϕ(i,r)  for i=1:Np]...)
ψᵣ = vcat([JacobiP(r,0.,0.,i) for i=0:Np-1]...)
println("|𝓥'*Iᵣ-ψᵣ| = ", norm(𝓥'*Iᵣ-ψᵣ))
