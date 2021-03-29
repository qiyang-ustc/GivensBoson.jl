

"""
    generate random particle-hole-symmetric positive-definite hamiltonian.
"""
function test_random_particle_hole_positive_definite_hamiltonian(N::Int)
    @assert mod(N,2)==0
    S = rand(N).+0.5
    q1,r1 = qr(rand(N,N))
    q2,r2 = qr(rand(N,N))
    S = diagm(vcat(S,S))
    q = [q1 q2;q2 q1]
    return transpose(q)*S*q
end

"""
    test if the transform V is canonical
"""
function test_canonicality(V::Matrix,ϵ::Float64)
    @assert size(V)[1]==size(V)[2]
    N = size(V)[1]÷2
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    @test norm(transpose(V)*η*V-η)<ϵ
end


"""
    test if the transform V is canonical with zeromodes
"""
function test_canonicality_zeromode(V::Matrix,ϵ::Float64)
    @assert size(V)[1]==size(V)[2]
    N = size(V)[1]÷2
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    test_η = transpose(V)*η*V
    for i = 1:2N
        if abs(test_η[i,i])<1E-10
            test_η[i,i] = i>N ? -1.0 : 1.0
        end
    end
    @test norm(test_η-η)<ϵ
end


"""
    test if the matrix is diagonalized
"""
function test_diag(s::Matrix,ϵ::Float64=1E-10)
    @test norm(s-diagm(diag(s)))<ϵ
end


"""
    test if S = transpose(V)*A*V
    this will hold if use normal and abnormal deal zeromodes
"""
function test_right(s::Matrix,V::Matrix,A::Matrix,ϵ::Float64)
    N = size(s)[1]÷2
    @test norm(transpose(V)*A*V-s)<ϵ
end

"""
    test if my_inverse(V^T)*S*my_inverse(V) = A
    this will not hold if use normal and abnormal deal zeromodes
"""
function test_full_right(s::Matrix,V::Matrix,A::Matrix,ϵ::Float64)
    N = size(s)[1]÷2
    @test norm(transpose(V)*A*V-s)<ϵ
    Vs = V[1:N,1:N]
    Us = V[1:N,1+N:2N]
    iV = [transpose(Vs) -transpose(Us);-transpose(Us) transpose(Vs)]
    @test norm(A-transpose(iV)*s*iV)<ϵ
end


"""
    OBC Schwinger Boson hamiltonian
"""

function schwinger_boson_hamiltonian(N::Int,bc::String)
    A = 1
    λ = 1
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
    return M
end