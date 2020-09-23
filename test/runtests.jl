using GalerkinEtena
using Test
@testset "Jacobi functions" begin
    x = [ -8.385864502265894e-01, -5.857254804559920e-01, -2.613290131006463e-01, 9.639069017068973e-02, 4.452559470863178e-01,  7.449277954410395e-01];
    w = [ 3.640915066793222e-02, 2.148950338971699e-01, 3.935534533149237e-01, 3.151164618818836e-01, 1.085195495207685e-01, 1.100558293610547e-02];
    @test JacobiGQ(3.14, 2., 5).points ≈ x atol=1e-12
    @test JacobiGQ(3.14, 2., 5).weights ≈ w atol=1e-12
    peval_check = [ 1.492369633161560e+01; 1.522362625780301e+03; 1.385512703514001e+04]
    @test JacobiP([1.,2,3],2.,3.14,5) ≈ peval_check atol=1e-11
    @test Legendre([0.65],8)[1][:,9] ≈ JacobiP([0.65], 0., 0., 8) atol=1e-12
end

@testset "integrate" begin
   q = JacobiGQ(0.,0.,10)
   @test integrate(x->x*x, q) ≈ 2/3. atol=1.e-12
end

@testset "assemble" begin
   m = Mesh1D([0., 0.5, 1.5, 3.0, 2.5], [[1; 2], [2; 3], [3; 5], [5; 4]])
   e2e = [1 2; 1 3; 2 4; 3 4 ]
   e2f = [1 1 ; 2 1; 2 1; 2 2 ]
   @test connect1D(m) == (e2e,e2f)
end

@testset "RK4" begin
   f=(x,y)->[x[1,1]  0. ; 0. 2*x[2,2]]
   sol=rk4(f,[1 2; 3 1.],0., 1.,1e-3)
   @test sol ≈ [exp(1.) 2 ; 3 exp(2.)] atol=1e-4
end
@testset "Discretization" begin
    P_check = [1, 10, 9, 19, 18, 28, 27, 37, 36, 46, 45, 55, 54, 64, 63, 73, 72, 82, 81, 90]
    M_check = [1, 9, 10, 18, 19, 27, 28, 36, 37, 45, 46, 54, 55, 63, 64, 72, 73, 81, 82, 90]
    m = Mesh1D(0., 1.,  10)
    ξ = JacobiGL(0., 0., 8)
    x, vmapM, vmapP = DGDiscretization(m, ξ)
    @test vmapM == M_check
    @test vmapP == P_check
end
@testset "Advec1D" begin
   u_final =[
-7.18408237e-11 1.98669331e-01 3.89418342e-01 5.64642473e-01 7.17356091e-01 8.41470985e-01 9.32039086e-01 9.85449730e-01 9.99573603e-01 9.73847631e-01 ;
 1.00240326e-02 2.08483569e-01 3.98631523e-01 5.72887296e-01 7.24303860e-01 8.46844716e-01 9.35624545e-01 9.87103975e-01 9.99230686e-01 9.71521222e-01 ;
 3.22757657e-02 2.30198224e-01 4.18943405e-01 5.90986635e-01 7.39469092e-01 8.58471251e-01 9.43248869e-01 9.90422132e-01 9.98110390e-01 9.66007136e-01 ;
 6.36452071e-02 2.60643087e-01 4.47249949e-01 6.16026367e-01 7.60243758e-01 8.74152629e-01 9.53211794e-01 9.94269412e-01 9.95688646e-01 9.57412917e-01 ;
 9.98334167e-02 2.95520207e-01 4.79425539e-01 6.44217687e-01 7.83326910e-01 8.91207360e-01 9.63558185e-01 9.97494987e-01 9.91664810e-01 9.46300088e-01 ;
 1.35890006e-01 3.30007713e-01 5.10969055e-01 6.71559672e-01 8.05377325e-01 9.07087125e-01 9.72634224e-01 9.99405466e-01 9.86333565e-01 9.33939658e-01 ;
 1.66933425e-01 3.59487507e-01 5.37709957e-01 6.94495607e-01 8.23593909e-01 9.19858121e-01 9.79450493e-01 9.99995264e-01 9.80673380e-01 9.22255142e-01 ;
 1.88835130e-01 3.80166032e-01 5.56340914e-01 7.10336239e-01 8.36012700e-01 9.28359973e-01 9.83696463e-01 9.99816080e-01 9.76076184e-01 9.13423211e-01 ;
 1.98669331e-01 3.89418342e-01 5.64642473e-01 7.17356091e-01 8.41470985e-01 9.32039086e-01 9.85449730e-01 9.99573603e-01 9.73847631e-01 9.09297427e-01 ;
  ]
  a = 2π
  α = 1.
  Np = 9
  K = 10
  ad = Advec1D(0., 2., K, Np)
  f = (x,t) -> rhs1D(ad, x, t, a, α)
  u = sin.(ad.x)
  tFinal = 10.
  xₘᵢₙ = minimum(abs.(ad.x[1, :]- ad.x[2,:]))
  CFL = 0.75
  dt = CFL ./ (2π) * xₘᵢₙ
  dt *= 0.5
  u = rk4(f, u, 0., tFinal, dt)
  using LinearAlgebra
  @test norm(u_final-u) < 1e-8
end

@testset "Maxwell1D" begin
E1=[ 
3.385869202032191e-05 4.526692168096674e-01  2.805527250091301e-01 -2.773789722196313e-01 -4.516163177950441e-01;
5.048915438902341e-02 4.660152878835446e-01  2.384202237558333e-01 -3.174143540595392e-01 -4.440420603967090e-01;
1.558310829364448e-01 4.767961262322121e-01  1.390023982779091e-01 -3.926887363404605e-01 -4.311861955849840e-01;
2.803310549782955e-01 4.532117708925854e-01 -1.500025472381084e-03 -4.565703654323375e-01 -4.125063513298208e-01;
3.799059249144263e-01 3.900235377047584e-01 -1.387689429760913e-01 -4.715240742407457e-01 -3.717755350612110e-01;
4.343660641963593e-01 3.193745821869441e-01 -2.350938134811906e-01 -4.593846076050846e-01 -3.173651740098322e-01;
4.527277511665667e-01 2.805202580137144e-01 -2.775039240661136e-01 -4.516945226029138e-01 -2.851350374345268e-01;
]
E2=[
-2.864352863680749e-01 -1.422059726180498e-02 -1.174430769418755e-01 -4.052750420841433e-03  1.264600143644388e-01;
-2.459592268901120e-01 -3.907942922338522e-02 -1.089944535364194e-01  9.414565077130191e-03  1.286277707624082e-01;
-1.417104103771918e-01 -6.968432878990674e-02 -9.104065224261977e-02  3.456205178427398e-02  1.237137786875264e-01;
 2.306455034341034e-03 -1.025914813939781e-01 -6.589925705058178e-02  6.403892268715325e-02  9.991861693725035e-02;
 5.132686808058452e-02 -1.274218633678368e-01 -4.092087963897927e-02  9.895379444297982e-02  5.844819670039177e-02;
 1.090834810744881e-02 -1.242356787886844e-01 -1.700048087478925e-02  1.208655509811961e-01  1.912069266436495e-02;
-1.567324312512593e-02 -1.170605142916132e-01 -4.210045083810123e-03  1.264286847805281e-01  2.925944370238357e-05;
]
E_check = [E1 E2]
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
  using LinearAlgebra
  @test norm(u[1:7,:]-E_check) < 1e-8
end
