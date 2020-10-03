include("./GivensBoson.jl")
using BenchmarkTools
using Main.GivensBoson 
# using Profile,ProfileView

N = 200
A = rand(2*N,2*N)
A = transpose(A)*A
origin_A = copy(A)
given_eigen_solver(A)
A = copy(origin_A)
given_eigen_solver(A)
# @profile given_eigen_solver(A)