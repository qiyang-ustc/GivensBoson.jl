using GivensBoson 
using LinearAlgebra

Nrange = 2:2:30
result = zeros(size(Nrange)[1],4)
for id in CartesianIndices(Nrange)
    N = Nrange[id]
    ϵ = 1E-11
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
    ϵ = 1E-4
    Lambda = diagm([λ for i = 1:N])
    M = [Lambda Q;transpose(Q) Lambda]
    A = M
    A = rand(2N,2N)
    A = transpose(A)*A
    origin_A = copy(A)
    S,V = recanonicalize(A,"AntiSymmetry")
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    result[id,1]=N
    result[id,2]=norm(transpose(V)*η*V-η)
    # result[id,3]=norm(S-diagm(S))
    result[id,4]=norm(transpose(V)*origin_A*V-diagm(S))
    # if N == 12
    #     vscodedisplay(V)
    # end
end
vscodedisplay(result)
