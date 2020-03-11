using Profile

 # This function tests the time performance
# and stability of Generalized Singular Value Decomposition (GSVD)
#
# First, we generate a m by n random matrix A and p by n random matrix B,
# Create a copy of A and B, respectively (original A and B will be overwritten)
#
# Then, we call gsvd function, return U, V, Q, C, S, R, k, l.
# Calculate and print time elapsed for computing GSVD.
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
function test(m, p, n)
    A = randn(m, n)
    B = randn(p, n)
    A_ = copy(A)
    B_ = copy(B)

    @time U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    e = eps(Float64)
    res_a = opnorm(U'*A_*Q - C*R, 1)/(max(m,n)*opnorm(A_, 1)*e)
    res_b = opnorm(V'*B_*Q - S*R, 1)/(max(m,n)*opnorm(B_, 1)*e)
    orthog_u = opnorm(I - U'*U, 1)/(m*e)
    orthog_v = opnorm(I - V'*V, 1)/(p*e)
    orthog_q = opnorm(I - Q'*Q, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return k, l
end

# Example from Julia v1.1.0 documentation
function test0()
    A = [1. 0.; 0. -1.]
    B = [0. 1.; 1. 0.]
    A_ = copy(A)
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    e = eps(Float64)
    res_a = opnorm(U'*A_*Q - C*R, 1)/(max(m,n)*opnorm(A_, 1)*e)
    res_b = opnorm(V'*B_*Q - S*R, 1)/(max(m,n)*opnorm(B_, 1)*e)
    orthog_u = opnorm(I - U'*U, 1)/(m*e)
    orthog_v = opnorm(I - V'*V, 1)/(p*e)
    orthog_q = opnorm(I - Q'*Q, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    U, V, Q, C, S, R, A_, B_
end

# Example 1 from Ch.4 of Jenny's Thesis
function test1()
    A  = Float64[1 2 1 0; 2 3 1 1; 3 4 1 2;4 5 1 3;5 6 1 4]
    A_ = copy(A)
    B = Float64[6 7 1 5;7 1 -6 13;-4 8 9 -2]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    @time U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    e = eps(Float64)
    res_a = opnorm(U'*A_*Q - C*R, 1)/(max(m,n)*opnorm(A_, 1)*e)
    res_b = opnorm(V'*B_*Q - S*R, 1)/(max(m,n)*opnorm(B_, 1)*e)
    orthog_u = opnorm(I - U'*U, 1)/(m*e)
    orthog_v = opnorm(I - V'*V, 1)/(p*e)
    orthog_q = opnorm(I - Q'*Q, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return U, V, Q, C, S, R, k, l
end

# Example 2 from Ch.4 of Jenny's Thesis
# It turns out this example is wrong
function test2()
    A = [1.0 2 1 0;2 3 1 1;3 4 1 2]
    A_ = copy(A)
    B = [4.0 5 1 3;5 6 1 4;6 7 1 5;7 1 -6 13]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    e = eps(Float64)
    res_a = opnorm(U'*A_*Q - C*R, 1)/(max(m,n)*opnorm(A_, 1)*e)
    res_b = opnorm(V'*B_*Q - S*R, 1)/(max(m,n)*opnorm(B_, 1)*e)
    orthog_u = opnorm(I - U'*U, 1)/(m*e)
    orthog_v = opnorm(I - V'*V, 1)/(p*e)
    orthog_q = opnorm(I - Q'*Q, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return U, V, Q, C, S, R, k, l
end

# Example 1 from MATLAB 2019b documentation
function test3()
    A = [1.0 6 11;2 7 12;3 6 13;4 9 14;5 10 15]
    A_ = copy(A)
    B = [8.0 1 6;3 5 7;4 9 2]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    e = eps(Float64)
    res_a = opnorm(U'*A_*Q - C*R, 1)/(max(m,n)*opnorm(A_, 1)*e)
    res_b = opnorm(V'*B_*Q - S*R, 1)/(max(m,n)*opnorm(B_, 1)*e)
    orthog_u = opnorm(I - U'*U, 1)/(m*e)
    orthog_v = opnorm(I - V'*V, 1)/(p*e)
    orthog_q = opnorm(I - Q'*Q, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    U, V, Q, C, S, R
end

# Example 2 from MATLAB 2019b documentation
function test4()
    A = [1.0 4 7 10 13;2 5 8 11 14;3 6 9 12 15]
    A_ = copy(A)
    B = [17.0 24 1 8 15;23 5 7 14 16;4 6 13 20 22;10 12 19 21 3;11 19 25 2 9]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    e = eps(Float64)
    res_a = opnorm(U'*A_*Q - C*R, 1)/(max(m,n)*opnorm(A_, 1)*e)
    res_b = opnorm(V'*B_*Q - S*R, 1)/(max(m,n)*opnorm(B_, 1)*e)
    orthog_u = opnorm(I - U'*U, 1)/(m*e)
    orthog_v = opnorm(I - V'*V, 1)/(p*e)
    orthog_q = opnorm(I - Q'*Q, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    U, V, Q, C, S, R
end

# Example from Ch. 93 LAPACK of Handbook of Linear Algebra, 2nd Edition
function test5()
    A = [1.0 2 3 1 5;0 3 2 0 2;1 0 2 1 0;0 2 3 0 -1;1 0 2 1 1;0 2 1 0 1]
    A_ = copy(A)
    B = [1.0 -2 2 1 1;0 3 0 0 0;1 -2 2 1 1;0 2 0 0 0;2 -4 4 2 2;1 3 2 1 1]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    e = eps(Float64)
    res_a = opnorm(U'*A_*Q - C*R, 1)/(max(m,n)*opnorm(A_, 1)*e)
    res_b = opnorm(V'*B_*Q - S*R, 1)/(max(m,n)*opnorm(B_, 1)*e)
    orthog_u = opnorm(I - U'*U, 1)/(m*e)
    orthog_v = opnorm(I - V'*V, 1)/(p*e)
    orthog_q = opnorm(I - Q'*Q, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    U, V, Q, C, S, R, k, l
end

# Example from Ch. 93 LAPACK of Handbook of Linear Algebra, 2nd Edition
# Same example from "A new preprocessing algorithm for the computation of
# the generalized singular value decomposition"
function test6()
    A = [1.0 2 3 1 5;0 3 2 0 2;1 0 2 1 0;0 2 3 0 -1;1 0 2 1 1;0 2 1 0 1]
    A_ = copy(A)
    B = [1.0 -2 2 1 1;0 3 0 0 0;1 -2 2 1 1;0 2 0 0 0;2 -4 4 2 2;1 3 2 1 1]
    B_ = copy(B)
    m,n = size(A)
    p = size(B)[1]
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    e = eps(Float64)
    res_a = opnorm(U'*A_*Q - C*R, 1)/(max(m,n)*opnorm(A_, 1)*e)
    res_b = opnorm(V'*B_*Q - S*R, 1)/(max(m,n)*opnorm(B_, 1)*e)
    orthog_u = opnorm(I - U'*U, 1)/(m*e)
    orthog_v = opnorm(I - V'*V, 1)/(p*e)
    orthog_q = opnorm(I - Q'*Q, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    U, V, Q, C, S, R, k, l
end

function test7()
    A = randn(3, 8)
    B = randn(4, 8)
    println("rank(B): ",rank(B))
    println("rank([A;B]): ", rank([A;B]))
    A_ = copy(A)
    B_ = copy(B)
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    return A_, B_, U, V, Q, C, S, R, k, l
end

function test8()
    A = randn(3, 8)
    B = randn(4, 8)
    println("rank(B): ",rank(B))
    println("rank([A;B]): ", rank([A;B]))
    A_ = copy(A)
    B_ = copy(B)
    U, V, Q, alpha, beta, R, k, l = gsvd(A, B, 0)
    return A_, B_, U, V, Q, alpha, beta, R, k, l
end

function test_full()
    m,p,n,r,rA,rB =5,4,3,3,1,2
    E = randn(r,n)
    A = randn(m1,rA)*randn(rA,r)*E
    B = randn(m2,rB)*randn(rB,r)*E
    show(A)
    show(B)
    println("rank of A: ", rank(A))
    println("rank of B: ", rank(B))
    println("rank of [A' B']': ", rank([A' B']'))
    U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
    show(R)
end

function test(dim)
    for i = 1:20
        m = dim*i*10
        n = m
        p = m
        A = randn(m, n)
        B = randn(p, n)
        A_ = deepcopy(A)
        B_ = deepcopy(B)

        @time U, V, Q, C, S, R, k, l = gsvd(A, B, 1)
        e = eps(Float64)
        res_a = opnorm(U'*A_*Q - C*R, 1)/(max(m,n)*opnorm(A_, 1)*e)
        res_b = opnorm(V'*B_*Q - S*R, 1)/(max(m,n)*opnorm(B_, 1)*e)
        orthog_u = opnorm(I - U'*U, 1)/(m*e)
        orthog_v = opnorm(I - V'*V, 1)/(p*e)
        orthog_q = opnorm(I - Q'*Q, 1)/(n*e)
        println("m: ", m)
        println("res_a: ", res_a)
        println("res_b: ", res_b)
        println("orthog_u: ", orthog_u)
        println("orthog_v: ", orthog_v)
        println("orthog_q: ", orthog_q)
        println("-----------------------------------")
    end
end

function test_profiling(m, p, n)
    maxI = 10
    t_pre = 0.0
    t_qr = 0.0
    t_csd = 0.0
    t_all = 0.0
    for i = 1:maxI
        A = randn(m, n)
        B = randn(p, n)
        ans = @timed gsvd(A, B, 1)
        t_pre = t_pre + ans[1][1]
        t_qr = t_qr + ans[1][2]
        t_csd = t_csd + ans[1][3]
        t_all = t_all + ans[2]
    end
    t_pre = t_pre/maxI
    t_qr = t_qr/maxI
    t_csd = t_csd/maxI
    t_all = t_all/maxI
    p_pre = 100*t_pre/t_all
    p_qr = 100*t_qr/t_all
    p_csd = 100*t_csd/t_all
    return t_pre, p_pre, t_qr, p_qr, t_csd, p_csd, t_all
end
