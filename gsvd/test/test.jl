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
#       res_a = ||U' * A_ * Q - C * R||/(max(m,n)* ||A_|| * eps)
#
#       res_b = ||V' * B_ * Q - S * R||/(max(p,n)* ||B_|| * eps)
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

    @time U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
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

function test0()
    A = [1. 0.; 0. -1.]
    B = [0. 1.; 1. 0.]
    A_ = copy(A)
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
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
    U, V, Q, C, S, R, A_, B_
end

function test1()
    A  = [1.0 2 1 0; 2 3 1 1; 3 4 1 2;4 5 1 3;5 6 1 4]
    A_ = copy(A)
    B = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
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

function test2()
    A = [1.0 2 1 0;2 3 1 1;3 4 1 2]
    A_ = copy(A)
    B = [4.0 5 1 3;5 6 1 4;6 7 1 5;7 1 -6 13]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
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

function test3()
    A = [1.0 2 1 0;-2 3 1 1;3 4 1 2]
    A_ = copy(A)
    B = [4.0 5 1 3;5 6 1 4;6 7 1 5;7 1 -6 13]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
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
    U, V, Q, C, S, R
end

function test4()
    A = [1.0 2 1 0;-2 3 1 1;3 4 1 2]
    A_ = copy(A)
    B = [4.0 5 1 3;5 -6 1 4;6 7 1 5;7 1 -6 13]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
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
    U, V, Q, C, S, R
end

function test5()
    A = [1.0 2 3 1 5;0 3 2 0 2;1 0 2 1 0;0 2 3 0 -1;1 0 2 1 1;0 2 1 0 1]
    A_ = copy(A)
    B = [1.0 -2 2 1 1;0 3 0 0 0;1 -2 2 1 1;0 2 0 0 0;2 -4 4 2 2;1 3 2 1 1]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, alpha, beta, R, k, l = gsvd(A, B, 0)
    # e = eps(Float64)
    # res_a = norm(U'*A_*Q - C*R, 1)/(max(m,n)*norm(A_, 1)*e)
    # res_b = norm(V'*B_*Q - S*R, 1)/(max(m,n)*norm(B_, 1)*e)
    # orthog_u = norm(I - U'*U, 1)/(m*e)
    # orthog_v = norm(I - V'*V, 1)/(p*e)
    # orthog_q = norm(I - Q'*Q, 1)/(n*e)
    # println("res_a: ", res_a)
    # println("res_b: ", res_b)
    # println("orthog_u: ", orthog_u)
    # println("orthog_v: ", orthog_v)
    # println("orthog_q: ", orthog_q)
    U, V, Q, alpha, beta, R, k, l
end

function test6()
    A = [1.0 2 3 1 5;0 3 2 0 2;1 0 2 1 0;0 2 3 0 -1;1 0 2 1 1;0 2 1 0 1]
    A_ = copy(A)
    B = [1.0 -2 2 1 1;0 3 0 0 0;1 -2 2 1 1;0 2 0 0 0;2 -4 4 2 2;1 3 2 1 1]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    # e = eps(Float64)
    # res_a = norm(U'*A_*Q - C*R, 1)/(max(m,n)*norm(A_, 1)*e)
    # res_b = norm(V'*B_*Q - S*R, 1)/(max(m,n)*norm(B_, 1)*e)
    # orthog_u = norm(I - U'*U, 1)/(m*e)
    # orthog_v = norm(I - V'*V, 1)/(p*e)
    # orthog_q = norm(I - Q'*Q, 1)/(n*e)
    # println("res_a: ", res_a)
    # println("res_b: ", res_b)
    # println("orthog_u: ", orthog_u)
    # println("orthog_v: ", orthog_v)
    # println("orthog_q: ", orthog_q)
    U, V, Q, C, S, R, k, l
end
