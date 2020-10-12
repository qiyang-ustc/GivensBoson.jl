# GivensBoson
## Diagonalize Quadratic Bosons system by using Givens Rotation.

GivensBoson.jl is a package can diagonize Quadratic Bosonic Hamiltonian with possible zero modes by using Givens Rotation.


## To Install:

```julia
julia>]add https://github.com/qiyang-ustc/GivensBoson.jl
```

## Note:
```
julia>  A = rand(2*N,2*N); 
        A = transpose(A)*A; # symmetric bosonic hamiltonian
        origin_A = copy(A)
        S,V = given_eigen_solver(A);
        η = diagm(vcat([1.0 for i=1:N],[-1.0 for i=1:N]));
```

Then see:
```
julia> @show norm(transpose(V)*η*V-η)
       @show norm(transpose(V)*origin_A*V-S)
```

Have fun~
