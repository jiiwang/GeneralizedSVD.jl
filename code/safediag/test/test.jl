include("../src/csd.jl")

using DataFrames
# using DelimitedFiles
using CSV

"""
and stability of CS Decomposition

First, we generate a m+p by m+p random matrix A
Then, we compute QR decomposition of A, A = QR
Afterwards, we choose Q1 = Q[1:m, 1:l]
Q2 = Q[m+1:m+p, 1:l], such that:

      Q1' * Q1 + Q2' * Q2 = I

Call csd function, return U, V, Z, C, S.
Calculate and print time elapsed for computing CSD.

Assume ||Q1' * Q1 + Q2' * Q2 - I|| ≈ eps,
compute:

      res_q1 = ||U * C * Z' - Q1||/(max(m,l)* ||Q1|| * eps)

      res_q2 = ||V * S * Z' - Q2||/(max(p,l)* ||Q2|| * eps)

      orthog_cs = ||C' * C + S' * S -I||/(max(m,p,l) * eps)

      orthog_u = ||U' * U - I||/(m * eps)

      orthog_v = ||V' * V - I||/(p * eps)

      orthog_z = ||Z' * Z - I||/(l * eps)

argument
m: num of rows in Q1
p: num of rows in Q2
l: num of cols in Q1, Q2

no return
"""
function test(m, p, l)
    if m + p < l
        println("input error, m + p >= l is required.")
        return
    else
        F = qr(randn(m+p, l))
        # Q1 = view(Q_, 1:m, 1:l)
        # Q2 = view(Q_, m+1:m+p, 1:l)

        Q1 = Matrix(F.Q)[1:m, :]
        Q2 = Matrix(F.Q)[m+1:m+p, :]

        e = eps(Float64)

        @time U, V, Z, C, S = csd(Q1, Q2, 1)

        res_q1 = opnorm(U * C * Z' - Q1, 1)/(max(m,l)*opnorm(Q1, 1)*e)
        res_q2 = opnorm(V * S * Z' - Q2, 1)/(max(p,l)*opnorm(Q2, 1)*e)
        orthog_cs = opnorm(C' * C + S' * S - I, 1)/(max(m,p,l)*e)
        orthog_u = opnorm(U' * U - I, 1)/(m*e)
        orthog_v = opnorm(V' * V - I, 1)/(p*e)
        orthog_z = opnorm(Z' * Z - I, 1)/(l*e)

        println("res_q1: ", res_q1)
        println("res_q2: ", res_q2)
        println("orthog_cs: ", orthog_cs)
        println("orthog_u: ", orthog_u)
        println("orthog_v: ", orthog_v)
        println("orthog_z: ", orthog_z)
        return res_q1, res_q2, orthog_cs, orthog_u, orthog_v, orthog_z
    end
end

function test()
    df = DataFrame(m = Int64[], p = Int64[], n = Int64[], res_q1 = Float64[],
    res_q2 = Float64[], orthog_cs = Float64[], orthog_u = Float64[],
    orthog_v = Float64[], orthog_z = Float64[])
    # case 1: m >=l and p >= l
    println("===============================================================")
    a1, b1, c1, d1, e1, f1 = test(60, 50, 40)
    push!(df, [60 50 40 a1 b1 c1 d1 e1 f1])
    println("---------------------------------------------------------------")
    a2, b2, c2, d2, e2, f2 = test(300, 250, 200)
    push!(df, [300 250 200 a2 b2 c2 d2 e2 f2])
    println("---------------------------------------------------------------")
    a3, b3, c3, d3, e3, f3 = test(900, 750, 600)
    push!(df, [900 750 600 a3 b3 c3 d3 e3 f3])
    println("---------------------------------------------------------------")
    a4, b4, c4, d4, e4, f4 = test(1500, 1250, 1000)
    push!(df, [1500 1250 1000 a4 b4 c4 d4 e4 f4])

    # case 2: m >=l and p < l
    println("===============================================================")
    a5, b5, c5, d5, e5, f5 = test(60, 40, 50)
    push!(df, [60 40 50 a5 b5 c5 d5 e5 f5])
    println("---------------------------------------------------------------")
    a6, b6, c6, d6, e6, f6 = test(300, 200, 250)
    push!(df, [300 200 250 a6 b6 c6 d6 e6 f6])
    println("---------------------------------------------------------------")
    a7, b7, c7, d7, e7, f7 =  test(900, 600, 750)
    push!(df, [900 600 750 a7 b7 c7 d7 e7 f7])
    println("---------------------------------------------------------------")
    a8, b8, c8, d8, e8, f8 = test(1500, 1000, 1250)
    push!(df, [1500 1000 1250 a8 b8 c8 d8 e8 f8])

    # case 3: m < l and p >=l
    println("===============================================================")
    a9, b9, c9, d9, e9, f9 = test(40, 60, 50)
    push!(df, [40 60 50 a9 b9 c9 d9 e9 f9])
    println("---------------------------------------------------------------")
    a10, b10, c10, d10, e10, f10 = test(200, 300, 250)
    push!(df, [200 300 250 a10 b10 c10 d10 e10 f10])
    println("---------------------------------------------------------------")
    a11, b11, c11, d11, e11, f11 = test(600, 900, 750)
    push!(df, [600 900 750 a11 b11 c11 d11 e11 f11])
    println("---------------------------------------------------------------")
    a12, b12, c12, d12, e12, f12 = test(1000, 1500, 1250)
    push!(df, [1000 1500 1250 a12 b12 c12 d12 e12 f12])

    # case 4: m < l and p < l
    println("===============================================================")
    a13, b13, c13, d13, e13, f13 = test(40, 50, 60)
    push!(df, [40 50 60 a13 b13 c13 d13 e13 f13])
    println("---------------------------------------------------------------")
    a14, b14, c14, d14, e14, f14 = test(200, 250, 300)
    push!(df, [200 250 300 a14 b14 c14 d14 e14 f14])
    println("---------------------------------------------------------------")
    a15, b15, c15, d15, e15, f15 = test(600, 750, 900)
    push!(df, [600 750 900 a15 b15 c15 d15 e15 f15])
    println("---------------------------------------------------------------")
    a16, b16, c16, d16, e16, f16 = test(1000, 1250, 1500)
    push!(df, [1000 1500 1250 a16 b16 c16 d16 e16 f16])

    # pwd() to check current directory
    # open("rand_matrix_test.csv", "w") do io
    #     writedlm(io, df)
    # end
    CSV.write("rand_matrix_test.csv", df, append=true)
    return df
end

"""
This function works similarly to the previous test function,
Instead of evaluating time performance and stability,
it returns the computed C, S.

argument
m: num of rows in Q1
p: num of rows in Q2
l: num of cols in Q1/Q2

returns C, S
"""
function testCS(m, p, l)
    if m + p < l
        println("input error, m + p >=l is required.")
        return
    else
        Q, R = qr(randn(m+p, m+p))
        Q_ = Q*Matrix(I,m+p,m+p)
        Q1 = view(Q_, 1:m, 1:l)
        Q2 = view(Q_, m+1:m+p, 1:l)
        U, V, Z, C, S = csd(Q1, Q2, 1)
        return C, S
    end
end

"""
This function works similarly to the previous test function,
Instead of evaluating time performance and stability,
it returns the computed U, V, Z, C, S.

argument
m: num of rows in Q1
p: num of rows in Q2
l: num of cols in Q1/Q2

returns U, V, Z, C, S
"""
function testProduct(m, p, l)
    if m + p < l
        println("input error, m + p >=l is required.")
        return
    else
        Q, R = qr(randn(m+p, m+p))
        Q_ = Q*Matrix(I,m+p,m+p)
        Q1 = view(Q_, 1:m, 1:l)
        Q2 = view(Q_, m+1:m+p, 1:l)
        return U, V, Z, C, S = csd(Q1, Q2, 1)
    end
end

# Example 1 in Chapter 3 of Jenny's thesis
function test1a()
    Q1 = [1/sqrt(7) 0 0 1/sqrt(3);
          1/sqrt(7) -2/sqrt(10) -1/2 -1/(2*sqrt(3));
          1/sqrt(7) -1/sqrt(10) 1/4  3/(4*sqrt(3));
          1/sqrt(7) 0 3/4 -3/(4*sqrt(3));
          1/sqrt(7) 0 0 0]
    Q2 = [1/sqrt(7) 1/sqrt(10) -1/4 -3/(4*sqrt(3));
          1/sqrt(7) 2/sqrt(10) -1/4 1/(4*sqrt(3))]
    U, V, Z, alpha, beta = csd(Q1, Q2, 0)
end

function test1b()
    Q1 = [1/sqrt(7) 0 0 1/sqrt(3);
          1/sqrt(7) -2/sqrt(10) -1/2 -1/(2*sqrt(3));
          1/sqrt(7) -1/sqrt(10) 1/4  3/(4*sqrt(3));
          1/sqrt(7) 0 3/4 -3/(4*sqrt(3));
          1/sqrt(7) 0 0 0]
    Q2 = [1/sqrt(7) 1/sqrt(10) -1/4 -3/(4*sqrt(3));
          1/sqrt(7) 2/sqrt(10) -1/4 1/(4*sqrt(3))]
    U, V, Z, C, S = csd(Q1, Q2, 1)
end

# Example 2 in Chapter 3 of Jenny's thesis
function test2()
    Q1 = [1/sqrt(7) 0 0 1/sqrt(3);
          1/sqrt(7) -2/sqrt(10) -1/2 -1/(2*sqrt(3));
          1/sqrt(7) -1/sqrt(10) 1/4  3/(4*sqrt(3))]
    Q2 = [1/sqrt(7) 0 3/4 -3/(4*sqrt(3));
          1/sqrt(7) 0 0 0;
          1/sqrt(7) 1/sqrt(10) -1/4 -3/(4*sqrt(3));
          1/sqrt(7) 2/sqrt(10) -1/4 1/(4*sqrt(3))]
    U, V, Z, C, S = csd(Q1, Q2, 1)
end
