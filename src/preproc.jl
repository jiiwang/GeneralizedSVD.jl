"""
    preproc(A, B) is the preprocessing step for computing the
    Generalized Singular Value Decomposition (GSVD), it computes:

```math
    A = U * RA * Q' ,    B = V * RB * Q'
```
where
- U is m-by-m orthogonal matrix,
- V is p-by-p orthogonal matrix,
- Q is n-by-n orthogonal matrix,
- RA and RB have the following structures:

(1) if m-k-l >= 0:

                   n-k-l   k    l                   n-k-l  k    l
       RA =     k ( 0     A12  A13 )     RB =   l ( 0      0   B13 )
                l ( 0      0   A23 )          p-l ( 0      0    0  )
            m-k-l ( 0      0    0  ),

(2) if m-k-l < 0:

                   n-k-l   k    l                   n-k-l  k    l
       RA =     k ( 0     A12  A13 )     RB =   l ( 0      0   B13 )
              m-k ( 0      0   A23 ),          p-l ( 0      0    0  )

where
- k-by-k matrix A12 and l-by-l matrix B13 are nonsingular upper triangular,
- A23 is l-by-l upper triangular if m-k-l >= 0,
  otherwise A23 is (m-k)-by-l upper trapezoidal,
- k+l is the effective numerical rank of the (m+p)-by-n matrix [A;B],
  l is the effective numerical rank of B.

# arguments:
A: m-by-n Float 32 or Float64 matrix,
B: p-by-n Float 32 or Float64 matrix.

# returns:
- U: m-by-m matrix,
- V: p-by-p matrix,
- Q: n-by-n matrix,
- k: interger, equals rank([A;B]) - rank(B),
- l: interger, equals rank(B),
- RA: m-by-n matrix, overwrites A in preproc!(),
- RB: p-by-n matrix, overwrites B in preproc!().
"""
function preproc!(A::StridedMatrix{T}, B::StridedMatrix{T}) where T
    m, n = size(A)
    p = size(B)[1]
    tola = max(m,n)*opnorm(A, 1)*eps(T)
    tolb = max(p,n)*opnorm(B, 1)*eps(T)

    # Step 1: QR decomposition with
    # col pivoting of B
    #               l   n-l
    # B * P = V * (B11  B12) l
    #             ( 0    0 ) p-l
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
    Q = F1.P

    # Determine the numerical rank of B
    l = 0
    for i = 1:min(p,n)
        if abs(B[i,i]) > tolb
            l += 1
        end
    end

    # Clean up B to make it upper triangular
    @views triu!(B[1:l,1:l])
    if p > l
        @views B[l+1:p,1:n] .= T(0.0)
    end

    # Step 2: RQ decomposition of (B11  B12)
    # (B11  B12) = (0  B22) * Z
    if p >= l && n != l
        B[1:l,:], tau1 = @views LAPACK.gerqf!(B[1:l,:])
        @views LAPACK.ormrq!('R', 'T', B[1:l,:], tau1, A)
        @views LAPACK.ormrq!('R', 'T', B[1:l,:], tau1, Q)
        # Clean up B
        @views B[1:l,1:n-l] .= T(0.0)
        @views triu!(B[1:l,n-l+1:n])
    end

    # Step 3: QR decomposition with
    # col pivoting of A11, if A11 is not empty
    #            n-l    l
    #       A = (A11  A12) m
    #
    #       A11 * P = U * (T11  T12)
    #                     ( 0    0 )
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
    @views triu!(A[1:k,1:k])
    if m > k
        @views A[k+1:m,1:n-l] .= T(0.0)
    end

    # Step 4: RQ decomposition of (T11  T12)
    # (T11  T12) = (0  T12) * Z
    if n-l > k
        A[1:k,1:n-l], tau2 = @views LAPACK.gerqf!(A[1:k,1:n-l])
        @views LAPACK.ormrq!('R', 'T', A[1:k,1:n-l], tau2, Q[1:n,1:n-l])
        # Clean up A
        @views A[1:k,1:n-l-k] .= T(0.0)
        @views triu!(A[1:k,n-l-k+1:n-l])
    end

    # Step 5: QR decomposition of A[k+1:m,n-l+1:n]
    if m > k
        # Weird, you have to add @views for submatrix input, otherwise,
        # it will not be overwritten.
        # A[k+1:m,n-l+1:n], tau3 = @views LAPACK.geqrf!(A[k+1:m,n-l+1:n])
        # @views LAPACK.ormqr!('R', 'N', A[k+1:m,n-l+1:n], tau3, U[:,k+1:m])
        @views F3 = qr!(A[k+1:m,n-l+1:n])
        if m - k > l
            U[:,k+1:m] = @views U[:,k+1:m]*(F3.Q*Matrix{T}(I,m-k,m-k))
        else
            U[:,k+1:m] = @views U[:,k+1:m]*Matrix{T}(F3.Q)
        end
        @views triu!(A[k+1:m,n-l+1:n])
    end

    return U, V, Q, k, l, A, B
end

# preproc() calls preproc!() before making a copy of input matrices A and B,
# thus A, B will not be overwritten
function preproc(A::StridedMatrix{T}, B::StridedMatrix{T}) where T
    return preproc!(copy(A),copy(B))
end
