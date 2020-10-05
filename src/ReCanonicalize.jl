using LinearAlgebra

inner_product(V1::Vector,V2::Vector,η::Matrix) = (transpose(V2)*η*V1)[1]
inner_product(V::Vector,η::Matrix) = (transpose(V)*η*V)[1]
function normalization!(V::Matrix,i::Int,η::Matrix)
    c = inner_product(V[:,i],η)
    if c < 1E-8
        nothing
    else
        V[:,i] .=  V[:,i]./sqrt(inner_product(V[:,i],η))
    end
end

function gram_schmit(V::Matrix,s::Int,e::Int,η::Matrix)
    # V_1' = V_1
    normalization!(V,s,η)
    for i = s+1:1:e
        for j = s:1:i-1
            V[:,i] .-= inner_product(V[:,j],V[:,i],η).*V[:,j]
        end
        normalization!(V,i,η)
    end
    return V
end

generate_metric(N::Int) = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))

function recanonicalize(h::Matrix,hamiltonian_type::String;tol=1E-10)
    N = Int(size(h)[1]/2)
    η = generate_metric(N)
    ih = η*h
    tS,tV = eigen(ih)
    @show norm(imag.(tS))
    @show norm(imag.(tV))
    S,V = copy(tS),copy(tV)
    for i = 1:N
        S[i] = tS[N+i]
        S[N+i] = tS[i]
        V[:,i] .= tV[:,N+i]
        V[:,N+i] .= tV[:,i]
    end #re sort eigen value from 0->large->0->small
    
    start_flag = 1
    end_flag = 1
    while end_flag <= N
        end_flag = start_flag
        for i = start_flag+1:N
            if abs(S[start_flag]-S[i]) < tol
                end_flag = i
            else
                break
            end
        end
        if start_flag!= end_flag
            V = gram_schmit(V,start_flag,end_flag,η)
        else
            normalization!(V,start_flag,η)
        end
        start_flag = end_flag+1
    end

    if hamiltonian_type == "Symmetry"
        V[1:N,N+1:2N] .= conj.(V[N+1:2N,1:N])
        V[N+1:2N,N+1:2N] .= conj.(V[1:N,1:N])
        S[N+1:2N] .= S[1:N]
    elseif hamiltonian_type == "AntiSymmetry"
        V[1:N,N+1:2N] .= conj.(-V[N+1:2N,1:N])
        V[N+1:2N,N+1:2N] .= conj.(V[1:N,1:N])
        S[N+1:2N] .= S[1:N]
    else
        error("We do not support this type of hamiltonian~~~~~hahahah")
    end
    return S,V
end

# N = 2
# η = generate_metric(N)
# h = [1 1 1 1;1 1 1 1; 1 1 1 1; 1 1 1 1]*1.0
# S,V, = recanonicalize(h)
# display(S)
# display(V)