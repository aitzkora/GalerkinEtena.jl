using GalerkinEtena
a = 2π
α = 1.
Np = 9
K = 10
ξ = RefGrid{1}(0., 2., Np)
ad = Advec{1}(0., 2., K, Np)
f = (x,t) -> rhs(ad, x, t, a, α)
u = sin.(ad.x)
tFinal = 10.
xₘᵢₙ = minimum(abs.(ad.x[1, :]- ad.x[2,:]))
CFL = 0.75
dt = CFL ./ (2π) * xₘᵢₙ
dt *= 0.5
u = rk4(f, u, 0., tFinal, dt)
