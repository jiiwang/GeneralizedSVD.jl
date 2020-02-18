using DelimitedFiles

A = readdlm("examples/full-rank/A2.csv", ',', Float64)
B = readdlm("examples/full-rank/B2.csv", ',', Float64)
A_ = copy(A);
B_ = copy(B);
U, V, Q, C, S, R, k, l = gsvd(A, B, 1);
