using DelimitedFiles

A = readdlm("a5.csv", ',', Float64)
B = readdlm("b5.csv", ',', Float64)

# U, V, Q, C, S, R, k, l = gsvd(A, B, 1);
