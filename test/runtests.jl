using GivensBoson
using Test
using LinearAlgebra
@testset "GivensBoson.jl" begin
    for N = 4:2:60
        ϵ = 1E-11
        η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
        A,B = rand(N,N),rand(N,N)
        A,B = map(x->transpose(x)*x,[A,B])
        A = [A B;transpose(B) A]
        A = A./maximum(A)/10
        r = [rand()+0.5 for i =1:N]
        r = vcat(r,r)
        r = diagm(r)
        A = r+A        
        origin_A = copy(A)
        S,V = given_eigen_solver(A,hamiltonian_type="Symmetry")
        η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
        @test norm(transpose(V)*η*V-η)<ϵ
        @test norm(S-diagm(diag(S)))<ϵ
        @test norm(transpose(V)*origin_A*V-S)<ϵ
    end
end

@testset "GivensBoson -Schwinger Boson PBC" begin
for N = 4:2:100
    A = 1
    λ = 1
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
    ϵ = 5*1E-10
    Lambda = diagm([λ for i = 1:N])
    M = [Lambda Q;transpose(Q) Lambda]
    A = M
    origin_A = copy(A)
    S,V = given_eigen_solver(A,zeromode=(N%4==0),hamiltonian_type="AntiSymmetry")
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    test_η = transpose(V)*η*V
    if N%4 == 0 && bc == "PBC"
        @test abs(sum(abs.(test_η-η))-4)< 1E-8 # 4 zero modes check
        for i = 1:2N
            if abs(test_η[i,i])<1E-10
                test_η[i,i] = i>N ? -1.0 : 1.0
            end
        end
    end
    @test norm(test_η-η)<ϵ
    @test norm(S-diagm(diag(S)))<ϵ
    @test norm(transpose(V)*origin_A*V-S)<ϵ
    end
end


@testset "GivensBoson -Schwinger Boson OBC" begin
for N = 4:2:256
    A = 1
    λ = 1
    bc = "OBC"
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
    ϵ = 5*1E-10
    Lambda = diagm([λ for i = 1:N])
    M = [Lambda Q;transpose(Q) Lambda]
    A = M
    origin_A = copy(A)
    S,V = given_eigen_solver(A,zeromode=(N%4==0),hamiltonian_type="AntiSymmetry")
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    test_η = transpose(V)*η*V
    @test abs(sum(abs.(test_η-η)))< 1E-8 # 4 zero modes check
    if N%4 == 0 && bc == "PBC"
        for i = 1:2N
            if abs(test_η[i,i])<1E-10
                test_η[i,i] = i>N ? -1.0 : 1.0
            end
        end
    end
    @test norm(test_η-η)<ϵ
    @test norm(S-diagm(diag(S)))<ϵ
    @test norm(transpose(V)*origin_A*V-S)<ϵ
    end
end


@testset "ReCanonicalize.jl" begin
    include("../src/ReCanonicalize.jl")

    A = 1
    λ = 1
    N = 20
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
    origin_A = copy(A)
    # S,V = given_eigen_solver(A)
    ϵ = 1E-6
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    origin_A = copy(A)
    S,V = recanonicalize(A,"AntiSymmetry")
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    temp = adjoint(V)*η*V-η
    for id in CartesianIndices(temp)
            if (abs(temp[id])-1)<ϵ
            temp[id] = 0
        end
    end
    # display(S)
    @test norm(temp)<ϵ
    # display(transpose(V)*origin_A*V-diagm(S))
    # display(V)
    # display(adjoint(V)*η*V-η)
    @test norm(adjoint(V)*origin_A*V-diagm(S))<ϵ
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