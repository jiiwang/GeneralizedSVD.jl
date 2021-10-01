using LinearAlgebra

preproc(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T = preproc!(copy(A),copy(B))

function preproc!(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    m, n = size(A)
    p = size(B)[1]
    tola = max(m,n)*opnorm(A, 1)*eps(T)
    tolb = max(p,n)*opnorm(B, 1)*eps(T)

    # Step 1: QR decomposition with
    # col pivoting of B
    #               l   n-l
    # B * P = V * (B11  B12) l
    #             ( 0    0 ) p-l
    #
    F1 = qr!(B, Val(true))
    if p > n
        V = F1.Q*Matrix{T}(I, p, p)
    else
        V = Matrix{T}(F1.Q)
    end

    # Step 2: Update A
    # A = A*P
    A = A*F1.P

    # Step 3: Set Q
    # Q = I_n * P
    Q = F1.P

    # Determine the numerical rank of B
    l = 0
    for i = 1:min(p,n)
        if abs(B[i,i]) > tolb
            l += 1
        end
    end

    # Clean up B to make it upper triangular
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

    # Step 2: RQ decomposition of (B11  B12)
    # (B11  B12) = (0  B22) * Z
    if p >= l && n != l
        # println("------RQ of (B11  B12)------")
        B[1:l,:], tau1 = @views LAPACK.gerqf!(B[1:l,:])
        @views LAPACK.ormrq!('R', 'T', B[1:l,:], tau1, A)
        @views LAPACK.ormrq!('R', 'T', B[1:l,:], tau1, Q)
        # Clean up B
        for i=1:l
            for j=1:n-l
                B[i, j] = 0.0
            end
        end
        for j = n-l+1:n
            for i = j-n+l+1:l
                B[i, j] = 0.0
            end
        end
    end

    # Step 3: QR decomposition with
    # col pivoting of A11, if A11 is not empty
    #            n-l    l
    #       A = (A11  A12) m
    #
    #       A11 * P = U * (T11  T12)
    #                     ( 0    0 )
    # println("------QR with col pivoting of A11------")
    F2 = @views qr!(A[:,1:n-l], Val(true))

    # Determine the numerical rank of A11
    k = 0
    for i = 1:min(m, n-l)
        if abs(A[i,i]) > tola
            k += 1
        end
    end


    if m > n-l
        U = F2.Q*Matrix{T}(I, m, m)
    else
        U = Matrix{T}(F2.Q)
    end

    A[1:m, n-l+1:n] = @views U' * A[1:m, n-l+1:n]

    Q[:,1:n-l] = @views Q[:,1:n-l]*F2.P

    # Clean up A to make it upper triangular
    for j = 1:k-1
        for i = j+1:k
            A[i,j] = 0.0
        end
    end

    if m > k
        for i = k+1:m
            for j = 1:n-l
                A[i,j] = 0.0
            end
        end
    end

    # Step 4: RQ decomposition of (T11  T12)
    # (T11  T12) = (0  T12) * Z
    if n-l > k
        # println("------RQ of (T11  T12)------")
        A[1:k,1:n-l], tau2 = @views LAPACK.gerqf!(A[1:k,1:n-l])
        @views LAPACK.ormrq!('R', 'T', A[1:k,1:n-l], tau2, Q[1:n,1:n-l])
        # Clean up A
        for i = 1:k
            for j = 1:n-l-k
                A[i,j] = 0.0
            end
        end

        for j = n-l-k+1:n-l
            for i = j-n+l+k+1:k
                A[i, j] = 0.0
            end
        end
    end

    # Step 5: QR decomposition of A[k+1:m,n-l+1:n]
    if m > k
        # Weird, you have to add @views for submatrix input, otherwise, it will not
        # be overwritten.
        # println("------QR of submatrix of A------")
        A[k+1:m,n-l+1:n], tau3 = @views LAPACK.geqrf!(A[k+1:m,n-l+1:n])
        @views LAPACK.ormqr!('R', 'N', A[k+1:m,n-l+1:n], tau3, U[:,k+1:m])

        for j = n-l+1:n
            for i = j-n+k+l+1:m
                A[i,j] = 0.0
            end
        end
    end

    return U, V, Q, k, l, A, B
end
