module GivensBoson
export given_eigen_solver
using LinearAlgebra
# macro output(printf)
#     if !isdefined(Main,:(runmode))
#         runmode = "Can Output"
#     end
#     if runmode == "Disable Output"
#         return nothing
#     else
#         return :($printf)
#     end
# end

function find_offdiagnoal_maximun(ih::Matrix;tol = 1E-13) #use 
    N = Int(size(ih)[1]/2)
    index = [1,1]
    temp_max = 0.0
    abnormal_index = [1,1]
    @inbounds for i = 1:2N
        @inbounds for j = i+1:2N
            t = abs(ih[i,j])
            if t>temp_max
                if abs(ih[i,j]/(ih[j,j]+ih[i,i]))>0.5
                    abnormal_index = [i,j]
                else
                    temp_max = t
                    index = [i,j]
                end    
            end
        end
    end
    if temp_max<tol 
        if abnormal_index!=[1,1]
            index = abnormal_index
            return ih[index...],index
        else
            return ih[1,2],[1,2]
        end
    end
    return ih[index...],index
end

function find_offdiagnoal_second_maximun(ih::Matrix) #use 
    N = Int(size(ih)[1]/2)
    index = [1,1]
    temp_max = 0.0
    temp_second = 0.0
    @inbounds for i = 1:2N
        @inbounds for j = i+1:2N
            t = abs(ih[i,j])
            if t>temp_max
                temp_max = t
                temp_second = temp_max
            end
        end
    end
    return temp_second
end

function given_transform!(H::Matrix,G::Matrix,i::Int,j::Int,temp_space::Matrix) # keep i<j
    N = Int(size(H)[1]/2)
    if i<=N && j>N
        given_abnormal_composite_rotation!(H,G,i,j,temp_space)
    else
        given_normal_composite_rotation!(H,G,i,j,temp_space)
    end
end

function deal_zeromodes!(G::Matrix;tol_max = 1000.0)
    N = Int(size(G)[1]/2)
    for i =1:2N
        if maximum(G[:,i])>tol_max
            G[:,i] .= sign.(G[:,i])*0.5
        end
    end
end

function given_eigen_solver(ih::Matrix ;max_iter=100000,tol=1E-13,zeromode=false)
    N = Int(size(ih)[1]/2)
    G = diagm([1.0 for i=1:2N])
    temp_space = zeros(Float64,2N,2N)
    @inbounds for iter = 1:max_iter
        temp_max,index = find_offdiagnoal_maximun(ih;tol)
        # print(temp_max,index,'\n')
        if abs(temp_max) < tol
            # @output print("iter<max_iter, converge!")
            if zeromode
                deal_zeromodes!(G)
            end
            return ih,G
        end
        i,j = index
        # @show temp_max,index
        given_transform!(ih,G,i,j,temp_space)
    end
    # @output print("Largest off-diagonal element is")
    print("Iteration doesnot converge!\n")
    if zeromode
        deal_zeromodes!(G)
    end
    # off_diag_ih = ih - diagm(diag(ih))
    # @output display(off_diag_ih)
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
    N = Int(size(H)[1]/2)
    t = -2H[i,j]/(H[j,j]+H[i,i])
    # if atanh(t)>1
    #     print(i,'\t',j,"\n")
    #     display(H)
    #     print("\n\n\n")
    # end
    if  abs(t) > 0.9999999
        # if abs(find_offdiagnoal_second_maximun(H))<1E-9 #zero mode
            H[i,j] = 0
            H[j,i] = 0
            H[i,i] = 0
            H[j,j] = 0
            H[i,:] .= 0 
            H[j,:] .= 0 
            H[:,i] .= 0 
            H[:,j] .= 0 
            G[:,i] .= sign.(G[:,i])*0.5
            G[:,j] .= sign.(G[:,i])*0.5
            return nothing
        # end
        t = sign(t)*0.9999999
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
    given_abnormal_rotation!(H,G,a+N,b,temp_space)
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

include("./ReCanonicalize.jl")
export recanonicalize

end