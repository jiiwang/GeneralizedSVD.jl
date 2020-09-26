using LinearAlgebra

function orthog_preproc(U, V, Q)
    println("---Testing the orthogonality in preprocessing step---")
    e = eps(Float64)
    dim_u = size(U)[2]
    dim_v = size(V)[2]
    dim_q = size(Q)[2]
    println("orthog_u: ", opnorm(U'*U-I,1)/(e*dim_u))
    println("orthog_v: ", opnorm(V'*V-I,1)/(e*dim_v))
    println("orthog_q: ", opnorm(Q'*Q-I,1)/(e*dim_q))
end

function orthog_qr(Q1, Q2)
    println("---Testing the orthogonality in qr step---")
    e = eps(Float64)
    dim_q = size(Q1)[2]
    Q = [Q1;Q2]
    println("orthog_q: ", opnorm(Q'*Q-I,1)/(e*dim_q))
end

function orthog_csd(U, V, Z)
    println("---Testing the orthogonality in csd step---")
    e = eps(Float64)
    dim_u = size(U)[2]
    dim_v = size(V)[2]
    dim_z = size(Z)[2]
    println("orthog_u: ", opnorm(U' * U - I, 1)/(dim_u*e))
    println("orthog_v: ", opnorm(V' * V - I, 1)/(dim_v*e))
    println("orthog_z: ", opnorm(Z' * Z - I, 1)/(dim_z*e))
end

function orthog_overall(U, V, Q)
    println("---Testing the overall orthogonality---")
    e = eps(Float64)
    dim_u = size(U)[2]
    dim_v = size(V)[2]
    dim_q = size(Q)[2]
    println("orthog_u: ", opnorm(U' * U - I, 1)/(dim_u*e))
    println("orthog_v: ", opnorm(V' * V - I, 1)/(dim_v*e))
    println("orthog_q: ", opnorm(Q' * Q - I, 1)/(dim_q*e))
end

function test_orthog()
    for i in 5:20:500
        println("-----------------------------------------")
        println("Input matrix dimension: m = n = p = ", i)
        A = randn(i, i)
        B = randn(i, i)
        gsvd(A, B, 1)
    end
end

function test_qr_orthog()
    # for i in 5:20:500
        i = 5000
        println("-----------------------------------------")
        println("Input matrix dimension: m = n = ", i)
        A = randn(i, i)
        F = qr(A)
        Q = Matrix{Float64}(F.Q)
        println(opnorm(Q'* Q - I, 1)/(size(Q, 2) * eps(Float64)))
    # end
end
