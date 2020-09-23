using GalerkinEtena
Np = 7
K = 10
pb = Maxwell1D(-2.0, 2.0, K, Np)
f = (x,t) -> rhs1D(pb, x, t)
E = sin.(π .* pb.x) .* (pb.x .< 0) 
H = zeros(Np, K)
u = [ E ; H ]
tFinal = 10.
xₘᵢₙ = minimum(abs.(pb.x[1, :]- pb.x[2,:]))
CFL = 1.
dt = CFL * xₘᵢₙ
u = rk4(f, u, 0., tFinal, dt)
