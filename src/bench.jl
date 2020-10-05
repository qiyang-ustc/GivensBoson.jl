using BenchmarkTools
using GivensBoson 

N = 400
A = rand(2*N,2*N)
A = transpose(A)*A
@btime given_eigen_solver(A)