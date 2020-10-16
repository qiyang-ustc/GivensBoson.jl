# using GivensBoson 
using LinearAlgebra
include("./GivensBoson.jl")
using Main.GivensBoson

Nrange = 4:40:300
result = zeros(size(Nrange)[1],4)
for id in CartesianIndices(Nrange)
    N = Nrange[id]
    ϵ = 1E-11
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
    ϵ = 1E-10
    Lambda = diagm([λ for i = 1:N])
    M = [Lambda Q;transpose(Q) Lambda]
    A = M
    origin_A = copy(A)
    # S,V = recanonicalize(A,"AntiSymmetry")
    S,V = given_eigen_solver(A,zeromode=true)

    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    result[id,1]=N
    if N%4 == 0
        t = transpose(V)*η*V
        print("Number of zero mode is: ",sum(abs.(t-η)),'\n')
        for i = 1:2N
            if abs(t[i,i])<1E-11
                t[i,i] = i>N ? -1.0 : 1.0
            end
        end
    end
    # display(t)
    result[id,2]=norm(t-η)
    result[id,3]=norm(S-diagm(diag(S)))
    result[id,4]=norm(transpose(V)*origin_A*V-S)
    # if N == 12
        # display(V)
        # display(transpose(V)*origin_A*V-S)
    # end
end
vscodedisplay(result)
# display(result)
