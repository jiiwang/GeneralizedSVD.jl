using LinearAlgebra

function preproc(A, B)
    m, n = size(A)
    p = size(B)[1]
    tola = max(m,n)*norm(A, 1)*eps(Float64)
    tolb = max(p,n)*norm(B, 1)*eps(Float64)

    # Step 1: QR decomposition with
    # col pivoting of B
    F1 = qr(B, Val(true))
    Q1 = Matrix(F1.Q)
    R1 = Matrix(F1.R)
    # formulate permutation matrix
    # F1.jpvt

    # compute effective rank of B
    # l = 0
    # for i = 1:min(p,n)
    #     if abs(R1[i,i]) >= tolb
    #         l += 1
    #     end
    # end
    l = rank(B)

    # Step 2: RQ decomposition of [B11  B12]
    R1_ = similar(R1, size(R1[1:l,:]))
    copyto!(R1_ , R1[1:l,:])
    # r1_m, r1_n = size(R1)
    r1_m = l
    r1_n = n
    # tau = similar(R1_, min(r1_m, r1_n))
    R1_, tau = LAPACK.gerqf!(R1_)
    Q = LAPACK.orgrq!(R1_, tau, length(tau))
    R_l = zeros(Float64, r1_m, r1_n-r1_m)
    R_r = UpperTriangular(R1_[1:r1_m,r1_n-r1_m+1:r1_n])
    R_b = [R_l R_r]
    Q_t = zeros(Float64, r1_n-r1_m, r1_n)
    Q_b = [Q_t; Q]
    # return R1[1:l,:], R_b, Q_b
    return R1[1:l,:], R_r, Q
    # Step 3: QR decomposition with
    # col pivoting of A11

    # Step 4: RQ decomposition of [A11 A12]

    # Step 5: QR decomposition of A23

end
