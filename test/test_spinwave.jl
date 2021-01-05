using GivensBoson,LinearAlgebra
begin
	export adj_mat,boundary_condition

    abstract type BoundaryType end
    abstract type Free <: BoundaryType end
    abstract type Periodic <: BoundaryType end
    abstract type Cylinder <: BoundaryType end
    struct BoundaryCondition{T<:BoundaryType} name::String end
    boundary_condition(name::String) = BoundaryCondition{eval(Meta.parse(name))}(name)


    """
    adj_mat is a adjacent matrix builder which will return a matrix which is simple to use in site-site coupling system.
    The result is a matrix A_{ij}
    A_{ij}=k means the j-th bond of site-i is point to site-k
    """
    function adj_mat(Lx::Int,Ly::Int,bc::BoundaryCondition{Periodic};D=2)
        N = Lx*Ly
        adj = zeros(Int,N,4) # up right down left
        n(i::Int,j::Int) = (i-1)*Ly+j
        function periodic_index(x::Int,Lx::Int) 
            if x == 0 
                return Lx
            elseif x < Lx+1 
                return x
            else
                return x-Lx
            end
        end
        pix = periodic_index

        # 
        for i = 1:Lx # row
            for j = 1:Ly # column
                id = n(i,j)
                adj[id,1] = n(i,pix(j-1,Ly)) 
                adj[id,2] = n(pix(i+1,Lx),j) 
                adj[id,3] = n(i,pix(j+1,Ly))
                adj[id,4] = n(pix(i-1,Lx),j)
            end
        end

        if D==1
            adj[:,2].=0
            adj[:,4].=0
        end
        return adj
    end

    function adj_mat(Lx::Int,Ly::Int,bc::BoundaryCondition{Free};D=2)
        adj = adj_mat(Lx,Ly,boundary_condition("Periodic"))
        n(i::Int,j::Int) = (i-1)*Ly+j
        for i = 1:Lx
            for j =1:Ly
                id = n(i,j)
                if j == 1
                    adj[id,1] = 0 
                end
                if j == Ly
                    adj[id,3] =0
                end

                if i == 1
                    adj[id,4] = 0 
                end
                if i == Lx
                    adj[id,2] =0
                end
            end
        end
        return adj
    end

    """ Lx = 5 ; Ly = 2
    |-1 2 3 4 5-|
    |-6 7 8 9 10-|
    -------
    """
    function adj_mat(Lx::Int,Ly::Int,bc::BoundaryCondition{Cylinder};D=2)
        adj = adj_mat(Lx,Ly,boundary_condition("Periodic"))
        n(i::Int,j::Int) = (i-1)*Ly+j
        for i = 1:Lx
            for j =1:Ly
                id = n(i,j)
                if j == 1
                    adj[id,1] = 0 
                end
                if j == Ly
                    adj[id,3] =0
                end
            end
        end
        return adj
    end

    function adj_mat(L::Int,bc::BoundaryCondition;D)
        if D == 2
            return adj_mat(L,L,bc,D=2)
        elseif D==1
            return adj_mat(L,1,bc,D=1)
        else
            error("Invalid D")
        end
    end
    # L = 6
    # bc = boundary_condition("Free")
    # bc = boundary_condition("Periodic")
    # @show adj_mat(L,bc)
end
function hamiltonian(adj::Matrix)
    N = size(adj)[1]
    h = zeros(2N,2N)
    for id = CartesianIndices(adj)
        if adj[id] != 0
            if id[1] > adj[id]
                h[id[1],id[1]] += 0.5
                h[adj[id],adj[id]] += 0.5
                h[id[1]+N,id[1]+N] += 0.5
                h[adj[id]+N,adj[id]+N] += 0.5

                h[id[1]+N,adj[id]] += 0.5
                h[adj[id]+N,id[1]] += 0.5
                h[id[1],adj[id]+N] += 0.5
                h[adj[id],id[1]+N] += 0.5
            end
        end
    end
    return h
end

function initialize_givens_eigen_solver(ih::Matrix;perturbation::Int = 1000,hamiltonian_type = "Symmetry")
    N = Int(size(ih)[1]/2)
    η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
    # S,V = eigen(η*ih)
    r = [i/N/1000 for i=1:N]
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

begin
	Lx = 4
	Ly = 3
	bc_type = "Cylinder"
    bc = boundary_condition(bc_type)
    N = Lx*Ly
end;

η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]))
h = hamiltonian(adj_mat(Lx,Ly,bc))
s,G = initialize_givens_eigen_solver(h,hamiltonian_type="Symmetry")
print("-------initialize_givens_eigen_solver--------------\n")
@show norm(transpose(G)*hamiltonian(adj_mat(Lx,Ly,bc))*G-s)
@show norm(transpose(G)^(-1)*s*G^(-1)-hamiltonian(adj_mat(Lx,Ly,bc)))
@show norm(transpose(G)*η*G-η)
print("----------------------\n")

s,G = given_eigen_solver(h,hamiltonian_type = "Symmetry")
vscodedisplay(s)
vscodedisplay(G)
print("-------givens_eigen_solver--------------\n")
@show norm(transpose(G)*hamiltonian(adj_mat(Lx,Ly,bc))*G-s)
@show norm(transpose(G)^(-1)*s*G^(-1)-hamiltonian(adj_mat(Lx,Ly,bc)))
@show norm(transpose(G)*η*G-η)
print("----------------------\n")