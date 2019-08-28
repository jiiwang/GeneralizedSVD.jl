function test1()
    E = randn(4, 5)
    # B is 5 by 5 and rank = 4
    B = randn(5, 4)*randn(4, 4)*E
    A = randn(3, 5)
    A, B
end

# test input A: rank = 3
# A = [2 4 1 0 0;0 3 0 0 0; 0 0 1 0 0.0]
# B = [3 2 4.0 6;1 0 5 7]
function test_rq(A)
    A_ = similar(A, size(A))
    copyto!(A_ , A)
    m, n = size(A)
    # tau = similar(A_, min(m, n))
    A_, tau = LAPACK.gerqf!(A_)
    R = Matrix(UpperTriangular(A_[1:m,n-m+1:n]))
    Q = LAPACK.orgrq!(A_, tau, length(tau))
    R_l = zeros(Float64, m, n-m)
    Q_t = zeros(Float64, n-m, n)
    return A, [R_l R], [Q_t;Q]
    # return A, R, Q
end

function test_tall()
    B = [3 2 1.0 0;6 4 2 0;9 6 3 0;1 1 0 1;1 1 0 1]
    A = [4.0 1 2 3;8 2 4 6;0 0 0 1;1 1 0 1]

    # A = [1.0 6 11;2 7 12;3 6 13;4 9 14;5 10 15]
    # B = [8.0 1 6;3 5 7;4 9 2;6 10 14]

    # A = [1.0 2 3 1 5;1 3 2 1 2;1 1 2 1 1;0 2 3 0 -1;1 0 2 1 1;0 2 1 0 1];
    # B = [1.0 -2 2 1 1;0 3 0 0 0;1 -2 2 1 1;0 2 0 0 0;2 -4 4 2 2;1 3 2 1 1];
    A_ = copy(A)
    B_ = copy(B)
    # return preproc(A_, B_)
    return U, V, Q, k, l, alpha, beta = preproc(A_, B_)
    # return piv, tau = preproc(A_, B_)
    # return A11, P11, Q11, R11 = preproc(A_, B_)
    # return A_, P1, Q2 = preproc(A_, B_)
end

function test_dggsvp3()
    B = [3 2 1.0 0;6 4 2 0;9 6 3 0;1 1 0 1;1 1 0 1]
    A = [4.0 1 2 3;8 2 4 6;0 0 0 1;1 1 0 1]

    # A = [1.0 6 11;2 7 12;3 6 13;4 9 14;5 10 15]
    # B = [8.0 1 6;3 5 7;4 9 2;6 10 14]

    # A = [1.0 2 3 1 5;1 3 2 1 2;1 1 2 1 1;0 2 3 0 -1;1 0 2 1 1;0 2 1 0 1];
    # B = [1.0 -2 2 1 1;0 3 0 0 0;1 -2 2 1 1;0 2 0 0 0;2 -4 4 2 2;1 3 2 1 1];
    A_ = copy(A)
    B_ = copy(B)
    return U, V, Q, k, l, alpha, beta = dggsvp3!(A_, B_)
end

function per()
    jpvt = [2 3 1]
    P2 = zeros(Float64, length(jpvt), length(jpvt))
    for i=1:length(jpvt)
        k = jpvt[i]
        P2[k,i] = 1.0
    end
    return P2
end
