include("../src/gsvd.jl")

# Example 3 in the special case note
function test3()
    A = [0.  1  2  2; 1  2  1  2; 1  1  1  2; 1  2  2  3];
    B = [0.  0  1  1; 1  2  2  3; 1  1  1  2; 1  1  1  2];
    n = size(A)[1]
    Q1, Q2, R, A23, B13, l = gsvd(A, B)
    println("Stability test on householderqr")
    e = eps(Float64)
    println("res_a23b13 = ", opnorm([Q1;Q2]*R - [A23; B13], 1)/(2*l*opnorm([A23; B13], 1)*e))
    println("orthog_q = ", opnorm(Q1'*Q1 + Q2'*Q2 -I, 1)/(l*e))
    return Q1, Q2, R
end

# Example 4 in the special case note
function test4()
    A =  [4.  2  2  0  0; 0  0  2  0  0; 4  4  4  0  2; 2  2  2  0  0; 2  4  4  0  2];
    B =  [4.  4  2  0  2; 2  2  0  0  2; 2  2  0  0  2; 4  2  2  0  0; 2  0  0  0  0];
    n = size(A)[1]
    Q1, Q2, R, A23, B13, l = gsvd(A, B)
    println("Stability test on householderqr")
    e = eps(Float64)
    println("res_a23b13 = ", opnorm([Q1;Q2]*R - [A23; B13], 1)/(2*l*opnorm([A23; B13], 1)*e))
    println("orthog_q = ", opnorm(Q1'*Q1 + Q2'*Q2 -I, 1)/(l*e))
    return Q1, Q2, R
end
