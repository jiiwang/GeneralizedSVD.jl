using LinearAlgebra

function preproc(A, B)
    m, n = size(A)
    p = size(B)[1]
    tola = max(m,n)*norm(A, 1)*eps(Float64)
    tolb = max(p,n)*norm(B, 1)*eps(Float64)

    # Step 1: QR decomposition with
    # col pivoting of B
    #            l  n-l
    # B*P = Q1*(B11 B12) l
    #          ( 0   0 ) p-l
    #
    F1 = qr!(B, Val(true))
    V = F1.Q*Matrix(1.0I, p, p)
    # R1 = F1.R
    P1 = F1.P
    # P1 = F1.P*Matrix(1.0I, n, n)

    A = A*P1
    # return A
    Q = Matrix(1.0I, n, n)*P1
    # Initialize U
    U = Matrix(1.0I, m, m)

    # Determine the numerical rank of B
    l = 0
    k = 0
    for i = 1:min(p,n)
        if abs(B[i,i]) > tolb
            l += 1
        end
    end

    for j = 1:l - 1
        for i = j + 1:l
            B[i, j] = 0.0
        end
    end

    if p > l
        for i = l+1:p
            for j = 1:n
                B[i, j] = 0.0
            end
        end
    end

    # Step 2: RQ decomposition of [B11  B12]
    # [B11 B12] = [0 B22] * Z
    #
    if p >= l && n != l
        B1 = similar(B, size(B[1:l,:]))
        copyto!(B1 , B[1:l,:])
        B[1:l,:], tau1 = LAPACK.gerqf!(B[1:l,:])
        R_r = Matrix(UpperTriangular(B[1:l,n-l+1:n]))

        # A = A * Z'
        LAPACK.ormrq!('R', 'T', B[1:l,:], tau1, A)
        # Q = Q * Z'
        LAPACK.ormrq!('R', 'T', B[1:l,:], tau1, Q)

        # Q_b = LAPACK.orgrq!(B[1:l,:], tau1, length(tau1))
        # Q_t = zeros(Float64, n-l, n)
        # Q2 = [Q_t; Q_b]
        B[1:l,1:n-l] = zeros(Float64, l, n-l)
        B[1:l,n-l+1:n] = R_r
        # return B1, B[1:l,:], Q2

        # A = A*Q2'
        # Q = Q*Q2'

        # A[1:l,1:l] = A[1:l,:]*Q_b'
        # Q[1:l,1:l] = Q[1:l,:]*Q_b'
    end

    # return B

    # Step 3: QR decomposition with
    # col pivoting of A11, if A11 is not empty
    #
    # A11 = similar(A, size(A[:,1:n-l]))
    # copyto!(A11, A[:,1:n-l])
    # F2 = qr!(A[:,1:n-l], Val(true))
    # P2 = F2.P
    # # P2 = F2.P*Matrix(1.0I,n-l,n-l)
    # U = U * (F2.Q*Matrix(1.0I, m, m))
    # R3 = F2.R

    # A[:,1:n-l], jpvt, tau1 = LAPACK.geqp3!(A[:,1:n-l])
    A[:,1:n-l], tau2, jpvt = LAPACK.geqp3!(A[:,1:n-l])
    # tau1 = convert(Array{Float64,1}, tau1)
    LAPACK.ormqr!('L', 'T', A[:,1:n-l], tau2, A[:,n-l+1:n])
    LAPACK.ormqr!('R', 'N', A[:,1:n-l], tau2, U)
    # Determine the numerical rank of A11
    for i = 1:min(m, n-l)
        if abs(A[i,i]) > tola
            k += 1
        end
    end

    println("k: ", k)
    # return A

    # A[1:n-l,1:n-l] = R3
    # A[n-l+1:m,1:n-l] = zeros(Float64, m-n+l, n-l)
    # A[:,n-l+1:n] = U' * A[:,n-l+1:n]
    P2 = zeros(Float64, length(jpvt), length(jpvt))
    for i=1:length(jpvt)
        k = jpvt[i]
        P2[k,i] = 1.0
    end
    # P3[1:n-l,1:n-l] = P2

    Q[1:n,1:n-l] = Q[1:n,1:n-l]*P2

    for j = 1:k - 1
        for i = j + 1:k
            A[i, j] = 0.0
        end
    end

    if m > k
        for i = m-k+1:m
            for j = 1:n-l
                A[i, j] = 0.0
            end
        end
    end
    # return A

    # return A11, jpvt, U, A[:,1:n-l]

    # return A11, P2, U, A[:,1:n-l]

    # Step 4: RQ decomposition of [A11 A12]
    if n-l > k
        # A1_ = similar(A, size(A[1:k,1:n-l]))
        # copyto!(A1_ , A[1:k,1:n-l])
        A[1:k,1:n-l], tau = LAPACK.gerqf!(A[1:k,1:n-l])
        R_r = Matrix(UpperTriangular(A[1:k,n-l-k+1:n-l]))
        # Q_b = LAPACK.orgrq!(A[1:k,1:n-l], tau, length(tau))
        # Q_t = zeros(Float64, n-l-k, n-l)
        # Q3 = [Q_t; Q_b]
        LAPACK.ormrq!('R', 'T', A[1:k,1:n-l], tau, Q[1:n,1:n-l])
        A[1:k,1:n-l-k] = zeros(Float64, k, n-l-k)
        A[1:k,n-l-k+1:n-l] = R_r
        # Q[1:n,1:n-l] = Q[1:n,1:n-l]*Q3'

    end

    # return A

    # Step 5: QR decomposition of A23
    if m > k
        # F3 = qr!(A[k+1:m,n-l+1:n])
        A[k+1:m,n-l+1:n], tau = LAPACK.geqrf!(A[k+1:m,n-l+1:n])
        U1 = similar(A, size(A[k+1:m,n-l+1:n]))
        copyto!(U1 , A[k+1:m,n-l+1:n])

        # U1 = F3.Q*Matrix(1.0I,m-k,m-k)
        LAPACK.orgqr!(U1, tau, length(tau))

        LAPACK.ormqr!('R', 'N', U1, tau, U[:,k+1:m])
        # U[:,k+1:m] = U[:,k+1:m]*U1

        for j = n-l+1:n
            for i = j-n+k+l+1:m
                A[i,j] = 0.0
            end
        end
    end

    return U, V, Q, k, l, A, B
end
