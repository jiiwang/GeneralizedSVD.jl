using Test
import GeneralizedSVD:GenSVD, gsvd
using LinearAlgebra
using DataFrames
using CSV

"""
    test(m, p, n)

This function tests the time performance
and backward stability of Generalized Singular Value Decomposition (GSVD)

First, we generate a m by n random matrix A and p by n random matrix B,
Create a copy of A and B, respectively (original A and B will be overwritten)

Then, we call gsvd function, return U, V, Q, C, S, R, k, l.
Calculate and print time elapsed for computing GSVD.

Compute:

      res_a = ||U' * A_ * Q - C * R||/(max(m,n)* ||A_|| * eps)

      res_b = ||V' * B_ * Q - S * R||/(max(p,n)* ||B_|| * eps)

      orthog_u = ||U' * U - I||/(m * eps)

      orthog_v = ||V' * V - I||/(p * eps)

      orthog_q = ||Q' * Q - I||/(n * eps)

# Arguments
- m: num of rows in A
- p: num of rows in B
- n: num of cols in A, B

no return
"""
function computemetric(F::GenSVD, A, B)
    e = eps(Float64)
    m, n = size(A)
    p = size(B)[1]
    res_a = opnorm(F.U'*A*F.Q - F.C*F.R, 1)/(max(m,n)*opnorm(A, 1)*e)
    res_b = opnorm(F.V'*B*F.Q - F.S*F.R, 1)/(max(p,n)*opnorm(B, 1)*e)
    orthog_cs = opnorm(I - F.C'*F.C - F.S'*F.S, 1)/((F.k+F.l)*e)
    orthog_u = opnorm(I - F.U'*F.U, 1)/(m*e)
    orthog_v = opnorm(I - F.V'*F.V, 1)/(p*e)
    orthog_q = opnorm(I - F.Q'*F.Q, 1)/(n*e)
    return res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q
end

function computemetric(F::GeneralizedSVD, A, B)
    e = eps(Float64)
    m, n = size(A)
    p = size(B)[1]
    res_a = opnorm(F.U'*A*F.Q - F.D1*F.R0, 1)/(max(m,n)*opnorm(A, 1)*e)
    res_b = opnorm(F.V'*B*F.Q - F.D2*F.R0, 1)/(max(p,n)*opnorm(B, 1)*e)
    orthog_cs = opnorm(I - F.D1'*F.D1 - F.D2'*F.D2, 1)/((F.k+F.l)*e)
    orthog_u = opnorm(I - F.U'*F.U, 1)/(m*e)
    orthog_v = opnorm(I - F.V'*F.V, 1)/(p*e)
    orthog_q = opnorm(I - F.Q'*F.Q, 1)/(n*e)
    return res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q
end

function test(m, p, n)
    A = randn(m, n)
    B = randn(p, n)
    F = gsvd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F, A, B)
    return F.k, F.l, res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q
end

function test()
    df = DataFrame(m = Int64[], p = Int64[], n = Int64[], k = Int64[], l = Int64[], res_a = Float64[],
    res_b = Float64[], orthog_cs = Float64[], orthog_u = Float64[],
    orthog_v = Float64[], orthog_q = Float64[])

    # case 1: m >= n & p >= n
    println("=======================Case 1============================")
    for i = 1:20
        a, b, c, d, e, f, g, h = test(60, 50, 40)
        push!(df, [60 50 40 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(300, 250, 200)
        push!(df, [300 250 200 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(900, 750, 600)
        push!(df, [900 750 600 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(1500, 1250, 1000)
        push!(df, [1500 1250 1000 a b c d e f g h])
        println("---------------------------------------------------------------")
    end

    # case 2: m >= n > p
    println("=======================Case 2================================")
    for i = 1:20
        a, b, c, d, e, f, g, h = test(60, 40, 50)
        push!(df, [60 40 50 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(300, 200, 250)
        push!(df, [300 200 250 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(900, 600, 750)
        push!(df, [900 600 750 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(1500, 1000, 1250)
        push!(df, [1500 1000 1250 a b c d e f g h])
        println("---------------------------------------------------------------")
    end

    # case 3: p >= n > m
    println("=======================Case 3==============================")
    for i = 1:20
        a, b, c, d, e, f, g, h = test(40, 60, 50)
        push!(df, [40 60 50 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(200, 300, 250)
        push!(df, [200 300 250 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(600, 900, 750)
        push!(df, [600 900 750 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(1000, 1500, 1250)
        push!(df, [1000 1500 1250 a b c d e f g h])
        println("---------------------------------------------------------------")
    end

    # case 4: n > m & n > p
    println("=====================Case 4============================")
    for i = 1:20
        a, b, c, d, e, f, g, h = test(20, 30, 60)
        push!(df, [20 30 60 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(200, 300, 600)
        push!(df, [200 300 600 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(400, 600, 1200)
        push!(df, [400 600 1200 a b c d e f g h])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g, h = test(1000, 1500, 3000)
        push!(df, [1000 1500 3000 a b c d e f g h])
        println("---------------------------------------------------------------")
    end

    CSV.write("gsvd_rand_matrix_test.csv", df, append=true)
    return df
end

# Example 1 from Ch.4 of Jenny's Thesis
function test1()
    A  = Float64[1 2 1 0; 2 3 1 1; 3 4 1 2;4 5 1 3;5 6 1 4]
    B = Float64[6 7 1 5;7 1 -6 13;-4 8 9 -2]
    @time F1 = gsvd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F1, A, B)
    println("--------Residual of gsvd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    println([A; B] ≈ [F1.U*F1.C; F1.V*F1.S]*F1.R*F1.Q')

    F2 = svd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F2, A, B)
    println("--------Residual of svd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return F1, F2
end

# Example 2 from Ch.4 of Jenny's Thesis
# It turns out this example is wrong

# This thus is used as Example 2 in working notes
function test2()
    A = [1. 2 1 0;2 3 1 1;3 4 1 2]
    B = [4. 5 1 3;5 6 1 4;6 7 1 5;7 1 -6 13]
    F1 = gsvd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F1, A, B)
    println("--------Residual of gsvd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    println([A; B] ≈ [F1.U*F1.C; F1.V*F1.S]*F1.R*F1.Q')

    F2 = svd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F2, A, B)
    println("--------Residual of svd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return F1, F2
end

function test2_wrapper()
    A = [1. 2 1 0;2 3 1 1;3 4 1 2];
    B = [4. 5 1 3;5 6 1 4;6 7 1 5;7 1 -6 13];
    return U, V, Q, D1, D2, k, l, R = gsvdwrapper(A, B)
end

# Example 1 in working notes
function test3()
    A = [1. 2 3 0; 5 4 2 1; 0 3 5 2; 2 1 3 3; 2 0 5 3];
    B = [1. 0 3 -1; -2 5 0 1; 4 2 -1 2];
    F1 = gsvd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F1, A, B)
    println("--------Residual of gsvd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    println([A; B] ≈ [F1.U*F1.C; F1.V*F1.S]*F1.R*F1.Q')

    F2 = svd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F2, A, B)
    println("--------Residual of svd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return F1, F2
end

function test3_wrapper()
    A = [1. 2 3 0; 5 4 2 1; 0 3 5 2; 2 1 3 3; 2 0 5 3];
    B = [1. 0 3 -1; -2 5 0 1; 4 2 -1 2];
    return U, V, Q, D1, D2, k, l, R = gsvdwrapper(A, B)
end

# Example 4 in working notes
function test4()
    A = [1. 4 2 3 0;3 4 0 -2 1;4 7 5 6 3];
    B = [1. 4 2 3 0;2 5 3 4 1;3 6 4 5 2;0 1 -1 3 1];
    F1 = gsvd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F1, A, B)
    println("--------Residual of gsvd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    println([A; B] ≈ [F1.U*F1.C; F1.V*F1.S]*F1.R*F1.Q')

    F2 = svd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F2, A, B)
    println("--------Residual of svd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return F1, F2
end

# Early exit after preprocessing
function test5()
    # rank(A) = 3, rank(B) = 1, rank([A;B]) = 4
    A = [1. 2 3 4;3 5 2 1;4 6 2 1];
    B = [0 0 3. 0;0 0 0 0;0 0 2. 0; 0 0 0 0];
    F1 = gsvd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F1, A, B)
    println("--------Residual of gsvd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    println([A; B] ≈ [F1.U*F1.C; F1.V*F1.S]*F1.R*F1.Q')

    F2 = svd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F2, A, B)
    println("--------Residual of svd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return F1, F2
end

# Early exit in CSD
function test6()
    # rank(A) = 4, rank(B) = 1, rank([A;B]) = 4
    A = [1.0 2 3 4;2 5 2 1;4 2 3 1;2 4 5 1];
    B = [1. 1 0 0;0 0 0 0; 0 0 0 0; 0 0 0 0];
    F1 = gsvd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F1, A, B)
    println("--------Residual of gsvd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    println([A; B] ≈ [F1.U*F1.C; F1.V*F1.S]*F1.R*F1.Q')

    F2 = svd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F2, A, B)
    println("--------Residual of svd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return F1, F2
end

# Example in the document GSVD.jl VS. gsvd.m
function test7()
    U = [1/2 0 sqrt(3)/2 0;
        0 -sqrt(3)/2 0 1/2;
        sqrt(3)/2 0 -1/2 0;
        0 1/2 0 sqrt(3)/2]
    V = [1/2 0 sqrt(3)/2 0;
        0 -sqrt(3)/2 0 1/2;
        sqrt(3)/2 0 -1/2 0;
        0 1/2 0 sqrt(3)/2]
    Q = [1/2 0 sqrt(3)/2 0;
        0 -sqrt(3)/2 0 1/2;
        sqrt(3)/2 0 -1/2 0;
        0 1/2 0 sqrt(3)/2]
    R = [0 1 0 0;0 0 1 0;0 0 0 1]
    C = [sqrt(3)/2 0 0;0 1/2 0;0 0 0;0 0 0]
    S = [1/2 0 0;0 sqrt(3)/2 0;0 0 1;0 0 0]
    A = U*C*R*Q'
    B = V*S*R*Q'
    F1 = gsvd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F1, A, B)
    println("--------Residual of gsvd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    println([A; B] ≈ [F1.U*F1.C; F1.V*F1.S]*F1.R*F1.Q')

    F2 = svd(A, B)
    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F2, A, B)
    println("--------Residual of svd---------")
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)
    return F1, F2
end

# function profiling_function(m, p, n)
#     maxI = 50
#     t_pre = 0.0
#     t_qr = 0.0
#     t_csd = 0.0
#     t_all = 0.0
#     for i = 1:maxI
#         println("---$i-th run---")
#         A = randn(m, n)
#         B = randn(p, n)
#         ans = @timed gsvd(A, B)
#         t_pre = t_pre + ans[1][1]
#         t_qr = t_qr + ans[1][2]
#         t_csd = t_csd + ans[1][3]
#         t_all = t_all + ans[2]
#     end
#     t_pre = t_pre/maxI
#     t_qr = t_qr/maxI
#     t_csd = t_csd/maxI
#     t_all = t_all/maxI
#     p_pre = 100*t_pre/t_all
#     p_qr = 100*t_qr/t_all
#     p_csd = 100*t_csd/t_all
#     return m, p, n, t_pre, p_pre, t_qr, p_qr, t_csd, p_csd, t_all
# end
#
# function test_profiling()
#     df = DataFrame(m = Int64[], p = Int64[], n = Int64[],
#     t_pre = Float64[], p_pre = Float64[],
#     t_qr = Float64[], p_qr = Float64[],
#     t_csd = Float64[], p_csd = Float64[],
#     t_all = Float64[])
#
#     # case 1: m >= n & p >= n
#     println("=======================Case 1============================")
#     t1a = profiling_function(1500, 1200, 1000)
#     push!(df, [t1a[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t1b = profiling_function(500, 500, 500)
#     push!(df, [t1b[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t1c = profiling_function(650, 310, 230)
#     push!(df, [t1c[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t1d = profiling_function(430, 610, 210)
#     push!(df, [t1d[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#
#     # case 2: m >= n > p
#     println("=======================Case 2================================")
#     t2a = profiling_function(1500, 1000, 1200)
#     push!(df, [t2a[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t2b = profiling_function(720, 220, 540)
#     push!(df, [t2b[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t2c = profiling_function(440, 180, 440)
#     push!(df, [t2c[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t2d = profiling_function(370, 290, 350)
#     push!(df, [t2d[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#
#     # case 3: p >= n > m
#     println("=======================Case 3==============================")
#     # Catch LAPACK Exception: DBDSDC did not converge
#     t3a = profiling_function(1000, 1500, 1200)
#     push!(df, [t3a[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t3b = profiling_function(250, 300, 300)
#     push!(df, [t3b[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t3c = profiling_function(360, 660, 600)
#     push!(df, [t3c[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t3d = profiling_function(130, 520, 480)
#     push!(df, [t3d[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#
#     # case 4: n > m & n > p
#     println("=====================Case 4============================")
#     # Catch LAPACK Exception: DBDSDC did not converge
#     # t4a = profiling_function(1000, 1200, 1500)
#     t4a = profiling_function(1300, 1200, 1500)
#     push!(df, [t4a[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     t4b = profiling_function(260, 600, 770)
#     push!(df, [t4b[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     # Might early exit
#     # t4c = profiling_function(370, 250, 700)
#     t4c = profiling_function(370, 500, 700)
#     push!(df, [t4c[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#     # Might early exit
#     # t4d = profiling_function(120, 120, 400)
#     t4d = profiling_function(220, 220, 400)
#     push!(df, [t4d[i] for i in 1:10])
#     println("---------------------------------------------------------------")
#
#     CSV.write("LAPACK-JuliaX.X-profiling.csv", df, append=true)
#     return df
# end
#
# function wrapper_profiling_function(m, p, n)
#     maxI = 20
#     t_dggsvp3 = 0.0
#     t_dtgsja = 0.0
#     t_all = 0.0
#     for i = 1:maxI
#         println("---$i-th run---")
#         A = randn(m, n)
#         B = randn(p, n)
#         ans = @timed gsvdwrapper(A, B)
#         t_dggsvp3 = t_dggsvp3 + ans[1][1]
#         t_dtgsja = t_dtgsja + ans[1][2]
#         t_all = t_all + ans[2]
#     end
#     t_dggsvp3 = t_dggsvp3/maxI
#     t_dtgsja = t_dtgsja/maxI
#     t_all = t_all/maxI
#     p_dggsvp3 = 100*t_dggsvp3/t_all
#     p_dtgsja = 100*t_dtgsja/t_all
#     return m, p, n, t_dggsvp3, p_dggsvp3, t_dtgsja, p_dtgsja, t_all
# end
#
# function test_wrapper_profiling()
#     df = DataFrame(m = Int64[], p = Int64[], n = Int64[],
#     t_dggsvp3 = Float64[], p_dggsvp3 = Float64[],
#     t_dtgsja = Float64[], p_dtgsja = Float64[],
#     t_all = Float64[])
#
#     # case 1: m >= n & p >= n
#     println("=======================Case 1============================")
#     t1a = wrapper_profiling_function(1500, 1200, 1000)
#     push!(df, [t1a[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t1b = wrapper_profiling_function(500, 500, 500)
#     push!(df, [t1b[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t1c = wrapper_profiling_function(650, 310, 230)
#     push!(df, [t1c[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t1d = wrapper_profiling_function(430, 610, 210)
#     push!(df, [t1d[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#
#     # case 2: m >= n > p
#     println("=======================Case 2================================")
#     t2a = wrapper_profiling_function(1500, 1000, 1200)
#     push!(df, [t2a[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t2b = wrapper_profiling_function(720, 220, 540)
#     push!(df, [t2b[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t2c = wrapper_profiling_function(440, 180, 440)
#     push!(df, [t2c[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t2d = wrapper_profiling_function(370, 290, 350)
#     push!(df, [t2d[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#
#     # case 3: p >= n > m
#     println("=======================Case 3==============================")
#     t3a = wrapper_profiling_function(1000, 1500, 1200)
#     push!(df, [t3a[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t3b = wrapper_profiling_function(250, 300, 300)
#     push!(df, [t3b[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t3c = wrapper_profiling_function(360, 660, 600)
#     push!(df, [t3c[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t3d = wrapper_profiling_function(130, 520, 480)
#     push!(df, [t3d[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#
#     # case 4: n > m & n > p
#     println("=====================Case 4============================")
#     t4a = wrapper_profiling_function(1300, 1200, 1500)
#     push!(df, [t4a[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t4b = wrapper_profiling_function(260, 600, 770)
#     push!(df, [t4b[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t4c = wrapper_profiling_function(370, 500, 700)
#     push!(df, [t4c[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#     t4d = wrapper_profiling_function(220, 220, 400)
#     push!(df, [t4d[i] for i in 1:8])
#     println("---------------------------------------------------------------")
#
#     CSV.write("Julia1.3-profiling.csv", df, append=true)
#     return df
# end
