module GivensBoson
export given_eigen_solver
using LinearAlgebra

"""
    find the maximum-abs elements which is off-diagonal of ih.
    [If abs(a)<tol. a can be seen as zero.]
    return the value and the row and column index. 
    and a. and If there is any off-diagonal element. 
    
    Input manual:
        if zeromode = true: this function will save zero mode and avoid to return imformation about them.
        if zeromode = false: we suppose there is no zero mode in this problem.
    
    Output Format:value, index,flag
            
        value: ih[index...]
        flag:
            flag = 1: there is a non-zero-mode element we need to report.
            flag = 2: there is a zero mode we need to report.
            flag = 3: there is nothing off-diagonal.
"""
function find_offdiagonal_maximun(ih::Matrix;tol = 1E-13,zeromode=false) #use 
    N = Int(size(ih)[1]/2)
    index = [1,1]
    temp_max = 0.0
    abnormal_index = [1,1]
    flag = 1 

    @inbounds for i = 1:2N
        @inbounds for j = i+1:2N # iterate up-triangle part.
            t = abs(ih[i,j])
            if t>temp_max
                if zeromode && abs(ih[i,j]/(ih[j,j]+ih[i,i]))>0.49
                    abnormal_index = [i,j]
                else
                    temp_max = t
                    index = [i,j]
                end
            end
        end
    end

    if temp_max < tol && abnormal_index==[1,1] # true: all non-zeromode elements are zero
        return 0.0,[1,1],3
    end

    if temp_max < tol && abnormal_index!=[1,1]
        index = abnormal_index
        return ih[index...],index,2
    end

    if temp_max > tol
        return ih[index...],index,1
    end 
end

function given_transform!(H::Matrix,G::Matrix,i::Int,j::Int,temp_space::Matrix) # keep i<j
    N = Int(size(H)[1]/2)
    if i<=N && j>N
        given_abnormal_composite_rotation!(H,G,i,j,temp_space)
    else
        given_normal_composite_rotation!(H,G,i,j,temp_space)
    end
end

"""
This function return a prepared G to make iteration steps converge much faster
"""
function initialize_givens_eigen_solver(ih::Matrix;perturbation::Int = 1000,hamiltonian_type = "Symmetry")
    N = Int(size(ih)[1]/2)
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    # S,V = eigen(η*ih)
    r = [i/N/perturbation for i=1:N]
    r = vcat(r,r)
    # th = V*diagm(S+sort(vcat(-r,r)))*V^(-1)  #add perturbation
    th = ih + diagm(r)

    S,G = eigen(η*th)
    for i = 1:2N  # normalization
        c = transpose(G[:,i])*η*G[:,i]
        G[:,i] .=  G[:,i]./sqrt(abs(c))
    end

    for i = 1:N  # d~d^{\dagger} construction
        G[:,i] .= G[:,i+N]  #Why this construnction raise too large error?
    end
    if hamiltonian_type == "Symmetry"
        G[1:N,N+1:2N].= G[N+1:2N,1:N]
        G[N+1:2N,N+1:2N].= G[1:N,1:N]
    else hamiltonian_type == "AntiSymmetry"
        G[1:N,N+1:2N].= -G[N+1:2N,1:N]
        G[N+1:2N,N+1:2N].= G[1:N,1:N]
    end
    # the previous steps would raise small error, we need to use given eigen solver to eliminate them.
    return transpose(G)*ih*G,G
end

export deal_zeromodes

"""
    A normal output result for given_eigen_solver is ih and V (may include zero mode).
    We have two mathod to deal with it.
    type=normal:    
        This function provide API to elimante zeromode by "math-well-defined rotation(physically-ill)"
   
            a = [1 -1; 1 1] and
            [0 0;0 4 ] = a^T*[1 -1;-1 1]*a
            [4 0;0 0 ] = a^T*[1  1; 1 1]*a

    
    type=abnormal:    
        This function provide API to elimante zeromode by "math-ill-defined rotation(physical-well)"
   
            a = [1 1; 1 1] and
            [0 0;0 0] = a^T*[1 -1;-1 1]*a
            a = [1 -1; -1 1]
            [0 0;0 0] = a^T*[1  1; 1 1]*a        
"""
deal_zeromodes(ih::Matrix,G::Matrix;tol=1E-13,type=:normal) = deal_zeromodes!(copy(ih),copy(G),tol=tol,type=type)

export deal_zeromodes!
function deal_zeromodes!(ih::Matrix,G::Matrix;tol=1E-13,type=:normal)
    if type == :normal
        return deal_zeromodes_normal!(ih,G,tol=tol)
    elseif type == :abnormal
        return deal_zeromodes_abnormal!(ih,G,tol=tol)
    else
        error("deal_zeromodes: Invalid type:$type")
    end
end


function deal_zeromodes_abnormal!(ih::Matrix,G::Matrix;tol=1E-13)
    for iter = 1:1000
        N = size(G)[1]
        temp_max,index,flag = find_offdiagonal_maximun(ih;tol=tol,zeromode=true)
        if flag == 3
            break
        end

        i,j = index
        sgn = sign(ih[i,i]/ih[i,j])

        tG = Matrix{Float64}(I,N,N)
        tG[i,j] = -1*sgn
        tG[j,i] = -1*sgn
        ih[i,j]=0
        ih[j,i]=0
        G .= G*tG

        if sgn<0
            ih[j,j]=0
            ih[i,i]=0
        else
            ih[i,i]=0
            ih[j,j]=0
        end
    end
    return ih,G
end

function deal_zeromodes_normal!(ih::Matrix,G::Matrix;tol=1E-13)
    for iter = 1:1000
        N = size(G)[1]
        temp_max,index,flag = find_offdiagonal_maximun(ih;tol=tol,zeromode=true)
        if flag == 3
            break
        end

        i,j = index
        sgn = sign(ih[i,i]/ih[i,j])

        tG = Matrix{Float64}(I,N,N)
        tG[i,j] = -1
        tG[j,i] = 1
        ih[i,j]=0
        ih[j,i]=0
        G .= G*tG

        if sgn<0
            ih[j,j]=4*ih[i,i]    
            ih[i,i]=0
        else
            ih[i,i]=4*ih[j,j]    
            ih[j,j]=0
        end
    end
    return ih,G
end


"""
    Givens Roation Eigen solver; Mathematical motivation and derivation can be found in notes dir:
"""
function given_eigen_solver(ih::Matrix ;max_iter=10000000,tol=1E-13,initialize=true,hamiltonian_type="AntiSymmetry",zeromode=false)
    N = Int(size(ih)[1]/2)
    if initialize
        ih, G = initialize_givens_eigen_solver(ih,hamiltonian_type=hamiltonian_type)
    else 
        G = diagm([1 for i =1:2N])
    end

    temp_space = zeros(Float64,2N,2N)
    @inbounds for iter = 1:max_iter
        temp_max,index,flag = find_offdiagonal_maximun(ih;tol=tol,zeromode=true)
        if abs(temp_max) < tol || flag!=1
            # if zeromode && flag == 2
            #     deal_zeromodes!(ih,G) # make ih diagonalized, then G will be not canonical.
            # end
            return ih,G
        end
        i,j = index
        given_transform!(ih,G,i,j,temp_space)
    end
    print("Iteration doesnot converge!\n")
    return ih,G
end

function given_normal_composite_rotation!(H::Matrix,G::Matrix,i::Int,j::Int,temp_space::Matrix)
    N = Int(size(H)[1]/2)
    given_normal_rotation!(H,G,i,j,temp_space)
    if i > N
        given_normal_rotation!(H,G,i-N,j-N,temp_space)
    else
        given_normal_rotation!(H,G,i+N,j+N,temp_space)
    end
end

function given_abnormal_rotation!(H::Matrix,G::Matrix,i::Int,j::Int,temp_space::Matrix)
    tol_sign(x) = abs(x)>1E-11 ? sign(x) : 0 
    N = Int(size(H)[1]/2)
    t = -2H[i,j]/(H[j,j]+H[i,i])    
    if  abs(t) > 0.9999999 || isnan(t)
            error("Givens Rotation fail!")
            return nothing
    end
    θ = 0.5*atanh(t)
    c,s = cosh(θ),sinh(θ)
    #@show c,s
    temp_G = temp_space
    # G .= G*temp_G
    @inbounds  for id = 1:2N
        temp_space[id,i] =  c*G[id,i] 
        temp_space[id,i] += s*G[id,j]
        temp_space[id,j] =  c*G[id,j] 
        temp_space[id,j] += s*G[id,i]
        G[id,i] = temp_space[id,i]
        G[id,j] = temp_space[id,j]
    end
    # H .= temp_G*H*temp_G
    # #effect of temp_G
    temp_G .= H # save memry
    @inbounds  for id = 1:2N
        temp_G[i,id] = c*H[i,id] 
        temp_G[i,id] += s*H[j,id]   
        temp_G[j,id] = c*H[j,id]
        temp_G[j,id] += s*H[i,id]
    end
    temp_G[:,i] .= temp_G[i,:]
    temp_G[:,j] .= temp_G[j,:]
    temp_G[i,i] =  c^2*H[i,i]+2*c*s*H[i,j]+s^2*H[j,j]
    temp_G[j,j] =  c^2*H[j,j]+2*c*s*H[i,j]+s^2*H[i,i]
    temp_G[i,j] = c*s*(H[i,i]+H[j,j])+(c^2+s^2)H[i,j]
    temp_G[j,i] = temp_G[i,j]
    H .= temp_G
end

function given_abnormal_composite_rotation!(H::Matrix,G::Matrix,i::Int,j::Int,temp_space::Matrix)
    N = Int(size(H)[1]/2)
    a,b = i,j-N
    given_abnormal_rotation!(H,G,a,b+N,temp_space)
   #given_abnormal_rotation!(H,G,a+N,b,temp_space)
end


function given_normal_rotation!(H::Matrix,G::Matrix,i::Int,j::Int,temp_space::Matrix;ϵ=1E-10)
    N = Int(size(H)[1]/2)
    if abs(H[i,i]-H[j,j])<ϵ
        θ = pi/4
    else
        t = 2H[i,j]/(H[j,j]-H[i,i])
        θ = 1/2*atan(t)
    end
    c = cos(θ)
    s = sin(θ)

    temp_G = temp_space
    
    # G .= G*temp_G
    @inbounds  for id = 1:2N
        temp_space[id,i] =  c*G[id,i] 
        temp_space[id,i] -= s*G[id,j]
        temp_space[id,j] =  c*G[id,j] 
        temp_space[id,j] += s*G[id,i]
        G[id,i] = temp_space[id,i]
        G[id,j] = temp_space[id,j]
    end
    # H .= transpose(temp_G)*H*temp_G
    # @show H[i,j]
    # effect of temp_G
    temp_G .= H # save memry
    @inbounds  for id = 1:2N
        temp_G[i,id] = c*H[i,id] 
        temp_G[i,id] -= s*H[j,id]   
        temp_G[j,id] = c*H[j,id]
        temp_G[j,id] += s*H[i,id]
    end
    temp_G[:,i] .= temp_G[i,:]
    temp_G[:,j] .= temp_G[j,:]
    temp_G[i,i] =  c^2*H[i,i]-2*c*s*H[i,j]+s^2*H[j,j]
    temp_G[j,j] =  c^2*H[j,j]+2*c*s*H[i,j]+s^2*H[i,i]
    temp_G[i,j] = 0.0
    temp_G[j,i] = 0.0
    # display(transpose(tt)*H*tt-temp_G)
    H .= temp_G
    # H .= transpose(tt)*H*tt
end

# ReCanonicalize is not polished!!!
# include("./ReCanonicalize.jl")
# export recanonicalize

end
