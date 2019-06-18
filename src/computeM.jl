using GalerkinEtena
JacobiP([1,2,3],2.,3.14,5)
ξ₅=JacobiGL(0.,0.,5)
q10=JacobiGQ(0.,0.,10)
𝓥 = hcat([JacobiP(ξ₅,0.,0.,i) for i=1:5]...)
𝓜  = inv(𝓥*𝓥')
x=q10.points
ϕₖ=lagrange(ξ)
ϕ=(i,x)->vander(x,6)'*ϕₖ[i,:]
res=[sum(ϕ(i,x).*ϕ(j,x).*q10.weights) for i=1:5 for j=1:5]
