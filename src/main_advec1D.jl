using GalerkinEtena
a = 2π
α = 1.
Np = 9
K = 10
ad = Advec1D(0., 2., Np, K)
f = (x,t) -> rhs1D(ad, x, t, a, α)
u = sin.(ad.x)
tFinal = 10.
xₘᵢₙ = minimum(abs.(ad.x[1, :]- ad.x[2,:]))
CFL = 0.75
dt = CFL ./ (2π) * xₘᵢₙ
dt *= 0.5
u = rk4(f, u, 0., tFinal, dt)
