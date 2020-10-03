using GivensBoson
using Test
using LinearAlgebra
@testset "GivensBoson.jl" begin
    N = 10
    ϵ = 1E-11
    A = rand(2*N,2*N)
    A = transpose(A)*A
    origin_A = copy(A)
    S,V = given_eigen_solver(A)
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    @test norm(transpose(V)*η*V-η)<ϵ
    @test norm(S-diagm(diag(S)))<ϵ
    @test norm(transpose(V)*origin_A*V-S)<ϵ
end
