# This function tests the time performance
# and stability of Quotient Singular Value Decomposition (QSVD)
#
# First, we generate a m by n random matrix A and p by n random matrix B,
# Create a copy of A and B, respectively (A and B will be overwritten)
#
# Then, we call qsvd function, return U, V, Q, C, S, R.
# Calculate and print time elapsed for computing QSVD.
#
# Compute:
#
#       res_a = ||U * A_ * Q' - C * R||/(max(m,n)* ||A_|| * eps)
#
#       res_b = ||V * B_ * Q' - S * R||/(max(p,n)* ||B_|| * eps)
#
#       orthog_u = ||U' * U - I||/(m * eps)
#
#       orthog_v = ||V' * V - I||/(p * eps)
#
#       orthog_q = ||Q' * Q - I||/(n * eps)
#
# argument
# m: num of rows in A
# p: num of rows in B
# n: num of cols in A, B
#
# no return
function test(m, n, p)
    A = randn(m, n)
    B = randn(p, n)
    A_ = copy(A)
    B_ = copy(B)

    @time  U, V, Q, C, S, R = qsvd(A, B)
    e = eps(Float64)
    res_a = norm(U'*A_*Q - C*R, 1)/(max(m,n)*norm(A_, 1)*e)
    res_b = norm(V'*B_*Q - S*R, 1)/(max(m,n)*norm(B_, 1)*e)
    orthog_u = norm(I - U'*U, 1)/(m*e)
    orthog_v = norm(I - V'*V, 1)/(p*e)
    orthog_q = norm(I - Q'*Q, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
end
