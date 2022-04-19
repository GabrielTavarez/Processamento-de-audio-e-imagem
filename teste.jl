A = Array{Any, 2}(undef ,9, 9);
A
A[1,1] = []
append!(A[1,1], 2)
fill!(A, [])
size(A[1,1])

