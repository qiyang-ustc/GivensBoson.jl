using GivensBoson
using Test
using Random
using LinearAlgebra
include("./testutil.jl")
include("./function_test.jl")

@testset "One Zero Mode Cases" begin
    A = [1 1;1 1]
    ϵ = 1E-12
    S,V = given_eigen_solver(copy(A),hamiltonian_type="Symmetry",zeromode=true)
    # do not deal zero modes. Canonical and Right
    test_canonicality(V,ϵ)
    test_right(S,V,A,ϵ)
    test_full_right(S,V,A,ϵ)

    deal_zeromodes!(S,V,type=:abnormal)

    # deal zero modes. Not canonical but diagonal.
    test_diag(S,ϵ)
    test_right(S,V,A,ϵ)
    #TODO: zeromode canonicality

    #TODO: if type=normal test
end

@testset "GivensBoson.jl - random case with partile-hole symmetry" begin
    for N = 4:2:40
        Random.seed!(42)
        ϵ = 1E-9
        @show N
        A = test_random_particle_hole_positive_definite_hamiltonian(N)
        origin_A = copy(A)
        S,V = given_eigen_solver(origin_A,hamiltonian_type="Symmetry")
        test_canonicality(V,ϵ)
        test_diag(S,ϵ)
        test_right(S,V,A,ϵ)
    end
end

@testset "GivensBoson -Schwinger Boson PBC- AntiSymmetry" begin
    for N = 4:2:60
        ϵ = 1E-9
        A = schwinger_boson_hamiltonian(N,"PBC")
        origin_A = copy(A)
        S,V = given_eigen_solver(A,hamiltonian_type="AntiSymmetry")

        test_canonicality_zeromode(V,ϵ)
        test_right(S,V,A,ϵ)
        test_diag(S,ϵ)
    end
end


@testset "GivensBoson -Schwinger Boson OBC- AntiSymmetry" begin
    for N = 4:2:64
        ϵ = 1E-9
        A = schwinger_boson_hamiltonian(N,"OBC")
        origin_A = copy(A)
        S,V = given_eigen_solver(A,hamiltonian_type="AntiSymmetry")

        test_canonicality_zeromode(V,ϵ)
        test_right(S,V,A,ϵ)
        test_diag(S,ϵ)
    end
end


# @testset "ReCanonicalize - GivensBoson Equal OBC test" begin
#     include("../src/ReCanonicalize.jl")
#     for N = 2:2:40
#         A = 1
#         λ = 1
#         bc = "OBC"
#         begin
#             Q = zeros(N,N)
#             for id in CartesianIndices((1:N,1:N))
#             if id[1]-id[2]==1
#                 Q[id] = 1.0*A
#             elseif id[2] - id[1] ==1
#                 Q[id] = -1.0*A
#             end
#             end

#             if bc == "PBC"
#             Q[N,1] = -1.0*A
#             Q[1,N] =  1.0*A
#             end
#             Q = Q/2.0
#             Q
#         end
#         Lambda = diagm([λ for i = 1:N])
#         M = [Lambda Q;transpose(Q) Lambda]
#         A = M
#         origin_A = copy(A)
#         # S,V = given_eigen_solver(A)
#         ϵ = 1E-8
    
#         S1,V1 = given_eigen_solver(origin_A)
#         S2,V2 = recanonicalize(A,"AntiSymmetry")
#         @test norm(S1-diagm(S2))<ϵ
#         @test norm(V1-V2)<ϵ
#     end
# end