# function test1()
#     E = randn(4, 5)
#     # B is 5 by 5 and rank = 4
#     B = randn(5, 4)*randn(4, 4)*E
#     A = randn(3, 5)
#     A, B
# end
#
# # test input A: rank = 3
# # A = [2 4 1 0 0;0 3 0 0 0; 0 0 1 0 0.0]
# # B = [3 2 4.0 6;1 0 5 7]
# function test_rq(A)
#     A_ = similar(A, size(A))
#     copyto!(A_ , A)
#     m, n = size(A)
#     # tau = similar(A_, min(m, n))
#     A_, tau = LAPACK.gerqf!(A_)
#     R = Matrix(UpperTriangular(A_[1:m,n-m+1:n]))
#     Q = LAPACK.orgrq!(A_, tau, length(tau))
#     R_l = zeros(Float64, m, n-m)
#     Q_t = zeros(Float64, n-m, n)
#     return A, [R_l R], [Q_t;Q]
#     # return A, R, Q
# end
#
# function test_tall()
#     B = [3 2 1.0 0;6 4 2 0;9 6 3 0;1 1 0 1;1 1 0 1]
#     A = [4.0 1 2 3;8 2 4 6;0 0 0 1;1 1 0 1]
#
#     # A = [1.0 6 11;2 7 12;3 6 13;4 9 14;5 10 15]
#     # B = [8.0 1 6;3 5 7;4 9 2;6 10 14]
#
#     # A = [1.0 2 3 1 5;1 3 2 1 2;1 1 2 1 1;0 2 3 0 -1;1 0 2 1 1;0 2 1 0 1];
#     # B = [1.0 -2 2 1 1;0 3 0 0 0;1 -2 2 1 1;0 2 0 0 0;2 -4 4 2 2;1 3 2 1 1];
#     A_ = copy(A)
#     B_ = copy(B)
#     # return preproc(A_, B_)
#     return U, V, Q, k, l, alpha, beta = preproc(A_, B_)
#     # return piv, tau = preproc(A_, B_)
#     # return A11, P11, Q11, R11 = preproc(A_, B_)
#     # return A_, P1, Q2 = preproc(A_, B_)
# end
#
# function test_dggsvp3()
#     B = [3 2 1.0 0;6 4 2 0;9 6 3 0;1 1 0 1;1 1 0 1]
#     A = [4.0 1 2 3;8 2 4 6;0 0 0 1;1 1 0 1]
#
#     # A = [1.0 6 11;2 7 12;3 6 13;4 9 14;5 10 15]
#     # B = [8.0 1 6;3 5 7;4 9 2;6 10 14]
#
#     # A = [1.0 2 3 1 5;1 3 2 1 2;1 1 2 1 1;0 2 3 0 -1;1 0 2 1 1;0 2 1 0 1];
#     # B = [1.0 -2 2 1 1;0 3 0 0 0;1 -2 2 1 1;0 2 0 0 0;2 -4 4 2 2;1 3 2 1 1];
#     A_ = copy(A)
#     B_ = copy(B)
#     return U, V, Q, k, l, alpha, beta = dggsvp3!(A_, B_)
# end
#
# function per()
#     jpvt = [2 3 1]
#     P2 = zeros(Float64, length(jpvt), length(jpvt))
#     for i=1:length(jpvt)
#         k = jpvt[i]
#         P2[k,i] = 1.0
#     end
#     return P2
# end
#
include("generatebyrank.jl")

function test1()
    A1 = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B1 = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]

    A2 = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B2 = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]

    A = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]
    m, n = size(A1)
    p = size(B1)[1]

    U1, V1, Q1, k1, l1, A1, B1 = @time dggsvp3!(A1, B1)
    println("k: ", k1)
    println("l: ", l1)
    U2, V2, Q2, k2, l2, A2, B2 = @time preproc(A2, B2)
    # A2, U, Q = preproc(A2, B2)

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

function test2()
    A = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]
    return preproc(A,B)
end

function test3(m, n, p)
    A = randn(m, n)
    B = randn(p, n)
    A1 = copy(A)
    B1 = copy(B)
    A2 = copy(A)
    B2 = copy(B)
    U1, V1, Q1, k1, l1, A1, B1 = @time dggsvp3!(A1, B1)
    println("k: ", k1)
    println("l: ", l1)
    U2, V2, Q2, k2, l2, A2, B2 = @time preproc(A2, B2)
    # A2, U, Q = preproc(A2, B2)

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


function preproctest(m, p, n)
    # params: row(A): row(B): col ratio

    maxI = 25
    maxJ = 5
    X = zeros(maxI)
    twrapper = zeros(maxI)
    tnative = zeros(maxI)

    for i = 1:maxI
        X[i] = 10*m*i
        time1 = 0.0
        time2 = 0.0
        for j = 1:maxJ
            A1 = randn(10*m*i, 10*n*i)
            B1 = randn(10*p*i, 10*n*i)
            A2 = deepcopy(A1)
            B2 = deepcopy(B1)
            # A_ = deepcopy(A1)
            # B_ = deepcopy(B1)
            ans1 = @timed dggsvp3!(A1, B1)
            ans2 = @timed preproc(A2, B2)
            time1 = time1 + ans1[2]
            time2 = time2 + ans2[2]
        end
            twrapper[i] = time1/maxJ
            tnative[i] = time2/maxJ

        # e = eps(Float64)
        # res_a1 = opnorm(ans1[1][1]'*A_*ans1[1][3] - ans1[1][4]*ans1[1][6], 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)
        # res_a2 = opnorm(ans2[1].U'*A_*ans2[1].Q - ans2[1].D1*ans2[1].R0, 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)
        println("m: ", X[i])
        # println("residual err of A by new gsvd: ", res_a1)
        # println("residual err of A by current gsvd in Julia: ", res_a2)
        println("--------------------------------------------------------")
    end

    plot(X, twrapper, color="red", linewidth=2.0, linestyle="-", label="wrapper preproc");
    plot(X, tnative, color="blue", linewidth=2.0, linestyle="--", label="native preproc");
    xlabel("# row of A");
    ylabel("elapsed time in seconds");
    title("cpu timing of preproc, m:p:n = $m:$p:$n");
    legend();
    show()
end


function test_profiling(m, n, p, k, l)


end

# test timing of preproc
# m = 300
# n = 300
# p = 200
# k = 80
# l = 200
function test_profiling1()
    A11 = generator(300, 100, 80);
    A12 = randn(300, 200);
    A = [A11 A12];

    println("rank of A11: ",rank(A11));

    F = svd(A11);


    tola = max(300,300)*opnorm(A, 1)*eps(Float64);
    tola1 = max(300,300)*eps(opnorm(A, 1));
    println(tola);
    println(tola1);

    println(diag(A)[1:100])

    k = 0;
    for i = 1:min(300, 100)
        if abs(A[i,i]) > tola1
            k += 1
        end
    end
    println("k: ", k)

    k1 = 0;
    for i = 1:size(F.S)
        if S[i] > tola
            k1 += 1
        end
    end

    println("k1: ", k+1)

    B = generator(200, 300, 200);
    println("l: ", rank(B));
    # @time preproc(A,B);
end
