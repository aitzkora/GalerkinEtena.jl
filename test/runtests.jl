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

@testset "integrate " begin
   q = JacobiGQ(0.,0.,10)
   @test integrate(x->x*x, q) ≈ 2/3. atol=1.e-12
end

@testset "assemble" begin
   m = Mesh1D([0., 0.5, 1.5, 3.0, 2.5], [[1; 2], [2; 3], [3; 5], [5; 4]])
   h = [1 2; 1 3; 2 4; 3 4 ]
   g = [1 1 ; 2 1; 2 1; 2 2 ]
   @test faces2Vertices(m) == (h,g)
end
