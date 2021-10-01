include("../src/csd.jl")

using DataFrames
# using DelimitedFiles
using CSV

"""
Timing and backward stability testing of the CS Decomposition (special case)


First, we generate an (m+l)-by-(m+l) random matrix A
Then, we compute QR decomposition of A, A = QR
Afterwards, we choose Q1 = Q[1:m, 1:l]
Q2 = Q[m+1:m+l, 1:l], such that:

      Q1' * Q1 + Q2' * Q2 = I

Call csd function, return U, V, Z, C, S.
Calculate and print time elapsed for computing CSD.

Assume ||Q1' * Q1 + Q2' * Q2 - I|| ≈ eps,
compute:

      res_q1 = ||U * C * Z' - Q1||/(l * ||Q1|| * eps)

      res_q2 = ||V * S * Z' - Q2||/(l * ||Q2|| * eps)

      orthog_cs = ||C' * C + S' * S -I||/(l * eps)

      orthog_u = ||U' * U - I||/(m * eps)

      orthog_v = ||V' * V - I||/(p * eps)

      orthog_z = ||Z' * Z - I||/(l * eps)

argument
m: num of rows of Q1
l: num of cols of Q1, Q2 and num of rows of Q2

no return
"""

function test(m, l)
    if m > l
        throw(ArgumentError("m <= l is required"))
    else
        F = qr(randn(m+l, l))

        Q1 = @views Matrix(F.Q)[1:m, :]
        Q2 = @views Matrix(F.Q)[m+1:m+l, :]

        e = eps(Float64)
        csd2by1(Q1, Q2)
        # @time U, V, Z, C, S, R1 = csd(Q1, Q2)
        # res_q1 = opnorm(U * C * Z' - Q1, 1)/(l*opnorm(Q1, 1)*e)
        # res_q2 = opnorm(V * S * Z' - Q2, 1)/(l*opnorm(Q2, 1)*e)
        # orthog_cs = opnorm(C' * C + S' * S - I, 1)/(l*e)
        # orthog_u = opnorm(U' * U - I, 1)/(m*e)
        # orthog_v = opnorm(V' * V - I, 1)/(l*e)
        # orthog_z = opnorm(Z' * Z - I, 1)/(l*e)
        # orthog_upperr22 = opnorm(R1, 1)/ (maximum(size(R1))*e)
        #
        # println("res_q1: ", res_q1)
        # println("res_q2: ", res_q2)
        # println("orthog_cs: ", orthog_cs)
        # println("orthog_u: ", orthog_u)
        # println("orthog_v: ", orthog_v)
        # println("orthog_z: ", orthog_z)
        # println("orthog_upperr22: ", orthog_upperr22)
        # return res_q1, res_q2, orthog_cs, orthog_u, orthog_v, orthog_z, orthog_upperr22
    end
end

function test()
    df = DataFrame(m = Int64[], l = Int64[], res_q1 = Float64[],
    res_q2 = Float64[], orthog_cs = Float64[], orthog_u = Float64[],
    orthog_v = Float64[], orthog_z = Float64[], orthog_upperr22 = Float64[])

    # case 1: m < l
    println("===============================================================")
    for i = 1:20
        a, b, c, d, e, f, g = test(40, 60)
        push!(df, [40 60 a b c d e f g])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g = test(200, 300)
        push!(df, [200 300 a b c d e f g])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g = test(600, 900)
        push!(df, [600 900 a b c d e f g])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g = test(1000, 1500)
        push!(df, [1000 1500 a b c d e f g])
        println("---------------------------------------------------------------")
    end

    # case 2: m = l
    println("===============================================================")
    for i = 1:20
        a, b, c, d, e, f, g = test(50, 50)
        push!(df, [50 50 a b c d e f g])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g = test(250, 250)
        push!(df, [250 250 a b c d e f g])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g = test(750, 750)
        push!(df, [750 750 a b c d e f g])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f, g = test(1250, 1250)
        push!(df, [1250 1250 a b c d e f g])
        println("---------------------------------------------------------------")
    end

    CSV.write("rand_matrix_test_no_refinement.csv", df, append=true)
    return df
end

"""
This function works similarly to the previous test function,
Instead of evaluating time performance and stability,
it returns the computed C, S.

argument
m: num of rows of Q1
l: num of cols of Q1, Q2 and num of rows of Q2

returns C, S
"""
function testCS(m, p, l)
    if m > l
        throw(ArgumentError("m <= l is required"))
    else
        F = qr(randn(m+l, m+l))
        Q = F.Q*Matrix(I,m+l,m+l)
        Q1 = view(Q, 1:m, 1:l)
        Q2 = view(Q, m+1:m+l, 1:l)
        U, V, Z, C, S = safeDiag(Q1, Q2)
        return C, S
    end
end

"""
This function works similarly to the previous test function,
Instead of evaluating time performance and stability,
it returns the computed U, V, Z, C, S.

argument
m: num of rows of Q1
l: num of cols of Q1, Q2 and num of rows of Q2

returns U, V, Z, C, S
"""
function testProduct(m, l)
    if m > l
        throw(ArgumentError("m <= l is required"))
    else
        F = qr(randn(m+l, m+l))
        Q = F.Q*Matrix(I,m+l,m+l)
        Q1 = view(Q, 1:m, 1:l)
        Q2 = view(Q, m+1:m+l, 1:l)
        return U, V, Z, C, S = safeDiag(Q1, Q2)
    end
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
    m = 3; l = 4;
    e = eps(Float64)

    U, V, Z, C, S = safeDiag(Q1, Q2)

    res_q1 = opnorm(U * C * Z' - Q1, 1)/(l*opnorm(Q1, 1)*e)
    res_q2 = opnorm(V * S * Z' - Q2, 1)/(l*opnorm(Q2, 1)*e)
    orthog_cs = opnorm(C' * C + S' * S - I, 1)/(l*e)
    orthog_u = opnorm(U' * U - I, 1)/(m*e)
    orthog_v = opnorm(V' * V - I, 1)/(l*e)
    orthog_z = opnorm(Z' * Z - I, 1)/(l*e)

    println("res_q1: ", res_q1)
    println("res_q2: ", res_q2)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_z: ", orthog_z)
    return U,V,Z,C,S
end

# Example 3 in the special case note
function test3()
    A = [0.  1  2  2; 1  2  1  2; 1  1  1  2; 1  2  2  3];
    B = [0.  0  1  1; 1  2  2  3; 1  1  1  2; 1  1  1  2];
    n = size(A)[1]
    Q1, Q2, R, A23, B13, l = gsvd(A, B)
    U, V, Z, C, S = safeDiag(Q1, Q2)

    e = eps(Float64)
    res_q1 = opnorm(U * C * Z' - Q1, 1)/(l*opnorm(Q1, 1)*e)
    res_q2 = opnorm(V * S * Z' - Q2, 1)/(l*opnorm(Q2, 1)*e)
    orthog_cs = opnorm(C' * C + S' * S - I, 1)/(l*e)
    orthog_u = opnorm(U' * U - I, 1)/(l*e)
    orthog_v = opnorm(V' * V - I, 1)/(l*e)
    orthog_z = opnorm(Z' * Z - I, 1)/(l*e)

    println("res_q1: ", res_q1)
    println("res_q2: ", res_q2)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_z: ", orthog_z)
    return U,V,Z,C,S
end

# Example 4 in the special case note
function test4()
    A =  [4.  2  2  0  0; 0  0  2  0  0; 4  4  4  0  2; 2  2  2  0  0; 2  4  4  0  2];
    B =  [4.  4  2  0  2; 2  2  0  0  2; 2  2  0  0  2; 4  2  2  0  0; 2  0  0  0  0];
    n = size(A)[1]
    Q1, Q2, R, A23, B13, l = gsvd(A, B)
    U, V, Z, C, S = safeDiag(Q1, Q2)

    e = eps(Float64)
    res_q1 = opnorm(U * C * Z' - Q1, 1)/(l*opnorm(Q1, 1)*e)
    res_q2 = opnorm(V * S * Z' - Q2, 1)/(l*opnorm(Q2, 1)*e)
    orthog_cs = opnorm(C' * C + S' * S - I, 1)/(l*e)
    orthog_u = opnorm(U' * U - I, 1)/(l*e)
    orthog_v = opnorm(V' * V - I, 1)/(l*e)
    orthog_z = opnorm(Z' * Z - I, 1)/(l*e)

    println("res_q1: ", res_q1)
    println("res_q2: ", res_q2)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_z: ", orthog_z)
    return U,V,Z,C,S
end
