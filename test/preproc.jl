include("wrapper/preprocWrapper.jl")
include("../src/preproc.jl")

# using Plots
# using DelimitedFiles
using DataFrames
using CSV

function test1()
    A1 = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B1 = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]

    A2 = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B2 = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]

    A = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]
    m, n = size(A1)
    p = size(B1)[1]

    U1, V1, Q1, k1, l1, A1, B1 = dggsvp3!(A1, B1)
    println("k: ", k1)
    println("l: ", l1)
    U2, V2, Q2, k2, l2, A2, B2 = preproc!(A2, B2)

    e = eps(Float64)
    println("Stability test on dggsvp3")
    println("res_a = ", opnorm(U1*A1*Q1' - A, 1)/(max(m,n)*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V1*B1*Q1' - B, 1)/(max(m,n)*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U1'*U1-I, 1)/(m*e))
    println("orthog_v = ", opnorm(V1'*V1-I, 1)/(p*e))
    println("orthog_q = ", opnorm(Q1'*Q1-I, 1)/(n*e))
    println("--------------------------------------")
    println("Stability test on preproc")
    println("res_a = ", opnorm(U2*A2*Q2' - A, 1)/(max(m,n)*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V2*B2*Q2' - B, 1)/(max(m,n)*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U2'*U2-I, 1)/(m*e))
    println("orthog_v = ", opnorm(V2'*V2-I, 1)/(p*e))
    println("orthog_q = ", opnorm(Q2'*Q2-I, 1)/(n*e))

end

# Example 3 in the special case note
function test2()
    A = [0.  1  2  2; 1  2  1  2; 1  1  1  2; 1  2  2  3];
    B = [0.  0  1  1; 1  2  2  3; 1  1  1  2; 1  1  1  2];
    n = size(A)[1]
    U, V, Q, k, l, A1, B1 = preproc(A, B)
    println("Stability test on preproc")
    e = eps(Float64)
    println("res_a = ", opnorm(U*A1*Q' - A, 1)/(n*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V*B1*Q' - B, 1)/(n*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U'*U-I, 1)/(n*e))
    println("orthog_v = ", opnorm(V'*V-I, 1)/(n*e))
    println("orthog_q = ", opnorm(Q'*Q-I, 1)/(n*e))
    return U, V, Q, k, l, A1, B1
end

# Example 4 in the special case note
function test3()
    A =  [4.  2  2  0  0; 0  0  2  0  0; 4  4  4  0  2; 2  2  2  0  0; 2  4  4  0  2];
    B =  [4.  4  2  0  2; 2  2  0  0  2; 2  2  0  0  2; 4  2  2  0  0; 2  0  0  0  0];
    n = size(A)[1]
    U, V, Q, k, l, A1, B1 = preproc(A, B)
    println("Stability test on preproc")
    e = eps(Float64)
    println("res_a = ", opnorm(U*A1*Q' - A, 1)/(n*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V*B1*Q' - B, 1)/(n*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U'*U-I, 1)/(n*e))
    println("orthog_v = ", opnorm(V'*V-I, 1)/(n*e))
    println("orthog_q = ", opnorm(Q'*Q-I, 1)/(n*e))
    return U, V, Q, k, l, A1, B1
end

function test(m, p, n)
    A = randn(m, n)
    B = randn(p, n)
    U, V, Q, k, l, A1, B1 = preproc(A, B)
    e = eps(Float64)
    res_a = opnorm(U*A1*Q' - A, 1)/(n*opnorm(A, 1)*e)
    res_b = opnorm(V*B1*Q' - B, 1)/(n*opnorm(B, 1)*e)
    orthog_u = opnorm(U'*U-I, 1)/(n*e)
    orthog_v = opnorm(V'*V-I, 1)/(n*e)
    orthog_q = opnorm(Q'*Q-I, 1)/(n*e)
    return res_a, res_b, orthog_u, orthog_v, orthog_q
end

function test()
    df = DataFrame(m = Int64[], p = Int64[], n = Int64[], res_a = Float64[],
    res_b = Float64[], orthog_u = Float64[],
    orthog_v = Float64[], orthog_q = Float64[])

    # case 1: m >= n & p >= n
    println("=======================Case 1============================")
    for i = 1:20
        a, b, c, d, e = test(60, 50, 40)
        push!(df, [60 50 40 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(300, 250, 200)
        push!(df, [300 250 200 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(900, 750, 600)
        push!(df, [900 750 600 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(1500, 1250, 1000)
        push!(df, [1500 1250 1000 a b c d e])
        println("---------------------------------------------------------------")
    end

    # case 2: m >= n > p
    println("=======================Case 2================================")
    for i = 1:20
        a, b, c, d, e = test(60, 40, 50)
        push!(df, [60 40 50 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(300, 200, 250)
        push!(df, [300 200 250 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(900, 600, 750)
        push!(df, [900 600 750 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(1500, 1000, 1250)
        push!(df, [1500 1000 1250 a b c d e])
        println("---------------------------------------------------------------")
    end

    # case 3: p >= n > m
    println("=======================Case 3==============================")
    for i = 1:20
        a, b, c, d, e = test(40, 60, 50)
        push!(df, [40 60 50 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(200, 300, 250)
        push!(df, [200 300 250 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(600, 900, 750)
        push!(df, [600 900 750 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(1000, 1500, 1250)
        push!(df, [1000 1500 1250 a b c d e])
        println("---------------------------------------------------------------")
    end

    # case 4: n > m & n > p
    println("=====================Case 4============================")
    for i = 1:20
        a, b, c, d, e = test(20, 30, 60)
        push!(df, [20 30 60 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(200, 300, 600)
        push!(df, [200 300 600 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(400, 600, 1200)
        push!(df, [400 600 1200 a b c d e])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e = test(1000, 1500, 3000)
        push!(df, [1000 1500 3000 a b c d e])
        println("---------------------------------------------------------------")
    end

    CSV.write("preproc_rand_matrix_test.csv", df, append=true)
    return df
end
