function test1()
    E = randn(4, 5)
    # B is 5 by 5 and rank = 4
    B = randn(5, 4)*randn(4, 4)*E
    A = randn(3, 5)
    A, B
end

function test_rq(A)
    A_ = similar(A, size(A))
    copyto!(A_ , A)
    m, n = size(A)
    tau = similar(A_, min(m, n))
    A_, tau = LAPACK.gerqf!(A_, tau)
    Q = LAPACK.orgrq!(A_, tau, length(tau))
    R = Matrix(UpperTriangular(A_[1:m,n-m+1:n]))
    R_l = zeros(Float64, m, n-m)
    Q_t = zeros(Float64, n-m, n)
    A, [R_l R], [Q_t;Q]
end
