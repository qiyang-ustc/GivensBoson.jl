# using GivensBoson
using Test
using LinearAlgebra
# @testset "GivensBoson.jl" begin
#     N = 10
#     ϵ = 1E-11
#     A = rand(2*N,2*N)
#     A = transpose(A)*A
#     origin_A = copy(A)
#     S,V = given_eigen_solver(A)
#     η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
#     @test norm(transpose(V)*η*V-η)<ϵ
#     @test norm(S-diagm(diag(S)))<ϵ
#     @test norm(transpose(V)*origin_A*V-S)<ϵ
# end

@testset "ReCanonicalize.jl" begin
    include("../src/ReCanonicalize.jl")

    A = 1
    λ = 1
    N = 10
    bc = "PBC"
    begin
        Q = zeros(N,N)
        for id in CartesianIndices((1:N,1:N))
            if id[1]-id[2]==1
                Q[id] = 1.0*A
            elseif id[2] - id[1] ==1
                Q[id] = -1.0*A
            end
        end
        
        if bc == "PBC"
            Q[N,1] = -1.0*A
            Q[1,N] =  1.0*A
        end
        Q = Q/2.0
        Q
    end
    Lambda = diagm([λ for i = 1:N])
    M = [Lambda Q;transpose(Q) Lambda]
    A = M

    ϵ = 1E-11
    origin_A = copy(A)
    S,V = recanonicalize(M,"AntiSymmetry")
    display(S)
    print("\n")
    display(V)
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    # display(transpose(V)*origin_A*V-diagm(S))
    @test norm(transpose(V)*η*V-η)<ϵ
    # display(transpose(V)*η*V-η)
    @test norm(transpose(V)*origin_A*V-S)<ϵ
end 