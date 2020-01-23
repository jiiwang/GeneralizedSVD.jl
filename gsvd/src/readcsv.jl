using DelimitedFiles

A = readdlm("examples/non-full-rank/A1.csv", ',', Float64)
B = readdlm("examples/non-full-rank/B1.csv", ',', Float64)

# U, V, Q, C, S, R, k, l = gsvd(A, B, 1);
