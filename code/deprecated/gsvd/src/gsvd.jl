# module Generalizedsvd

using LinearAlgebra
# using Main.CSD

import Base.show
include("preproc.jl")
include("../../safediag/src/safeDiag.jl")

struct GSVD{T,S} <: Factorization{T}
    U::S
    V::S
    Q::S
    D1::S
    D2::S
    k::Int
    l::Int
    R0::S
    function GSVD{T,S}(U::AbstractMatrix{T}, V::AbstractMatrix{T}, Q::AbstractMatrix{T},
                                 D1::AbstractMatrix{T}, D2::AbstractMatrix{T}, k::Int, l::Int, R0::AbstractMatrix{T}) where {T,S}
        new(U, V, Q, D1, D2, k, l, R0)
    end

    function show(io::IO, mime::MIME{Symbol("text/plain")}, F::GSVD{<:Any,<:AbstractArray})
        summary(io, F); println(io)
        println(io, "U factor:")
        show(io, mime, F.U)
        println(io, "\nV factor:")
        show(io, mime, F.V)
        println(io, "\nQ factor:")
        show(io, mime, F.Q)
        println(io, "\nD1 factor:")
        show(io, mime, F.D1)
        println(io, "\nD2 factor:")
        show(io, mime, F.D2)
        println(io, "\nGeneralized singular values:")
        show(io, mime, sqrt.(diag(F.D1'*F.D1) ./ diag(F.D2'*F.D2)))
        println(io, "\nk factor:")
        show(io, mime, F.k)
        println(io, "\nl factor:")
        show(io, mime, F.l)
        println(io, "\neffective numerical rank of the matrix [A; B]:")
        show(io, mime, F.k + F.l)
        println(io, "\nR0 factor:")
        show(io, mime, F.R0)
    end
end

function GSVD(U::AbstractMatrix{T}, V::AbstractMatrix{T}, Q::AbstractMatrix{T},
                        D1::AbstractMatrix{T}, D2::AbstractMatrix{T}, k::Int, l::Int, R0::AbstractMatrix{T}) where T
    GSVD{T,typeof(U)}(U, V, Q, D1, D2, k, l, R0)
end


"""
This function computes the generalized singular value
decomposition (GSVD) of an m-by-n matrix A and p-by-n

```math
     A = U * D1 * R0 * Q',    B = V * D2 * R0 * Q'
```

 where U, V and Q are orthogonal matrices. Let k+l be the effective numerical rank of the matrix (A',B')',
 then R is a (k+l)-by-n matrix of structure [0 R] where R is (k+l)-by-(k+l) and is nonsingular upper triangular matrix,
 D1 and D2 are m-by-(k+l) and p-by-(k+l) "diagonal" matrices and
 have one of the following structures:

 If m >= k+l,

        D1 =      k   l
            k  (  I   0  )
            l  (  0   C  )
        m-k-l  (  0   0  )

        D2 =      k   l
            l  (  0   S  )
          p-l  (  0   0  )


 where

      C = diag( alpha( k+1 ), ... , alpha( k+l ) ),
      S = diag( beta( k+1 ), ..., beta( k+l ) ),
      C**2 + S**2 = I


 (2) If m < k+l,

       D1 =         k   m-k  k+l-m
                k ( I    0     0  )
              m-k ( 0    C     0  )

        D2 = 	      k  m-k   k+l-m
              m-k ( 0   S     0  )
            k+l-m ( 0   0     I  )
            p-k-l ( 0   0     0  )

 where

      C = diag( alpha(k+1), ... , alpha(m) ),
	    S = diag( beta(k+1), ..., beta(m) ),
      C**2 + S**2 = I


arguments:
A: m by n matrix
B: p by n matrix
# option: 0 or 1

returns U, V, Z, alpha, beta, R0, k, l
# if option = 0
# U, V, Z, D1, D2, R, k, l if option = 1
"""
gsvd(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T = gsvd!(copy(A),copy(B))

function gsvd!(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    m, n = size(A)
    p = size(B)[1]

    # Step 1:
    # preprocess A, B
    U, V, Q, k, l, A, B = preproc!(A, B)
    # profiling
    # ans1 = @timed preproc!(A, B)
    # U = ans1[1][1]
    # V = ans1[1][2]
    # Q = ans1[1][3]
    # k = ans1[1][4]
    # l = ans1[1][5]
    # A = ans1[1][6]
    # B = ans1[1][7]

    # Step 2:
    # QR decomposition of A23 and B13
    # A23 is upper triangular

    if m-k-l >= 0
        @views A23 = A[k+1:k+l, n-l+1:n]
    # A23 is upper trapezoidal
    else
        @views A23 = A[k+1:m, n-l+1:n]
    end
    @views B13 = B[1:l, n-l+1:n]

    # Quick exit when A23 doesn't exist
    if size(A23)[1] <= 0
        @views R0 = [A[1:k,:];B[1:l,:]]
        @views CS = [Matrix{T}(I, k+l, k+l);zeros(T, m+p-k-l, k+l)]
        @views C = CS[1:m, :]
        @views S = CS[m+1:m+p, :]
        return GSVD(U, V, Q, C, S, k, l, R0)
    end

    Q1, Q2, R23 = householderqr(A23, B13)
    # profiling
    # ans2 = @timed householderqr(A23, B13)
    # Q1 = ans2[1][1]
    # Q2 = ans2[1][2]
    # R23 = ans2[1][3]

    # Step 3:
    # CSD of Q1, Q2
    U1, V1, Z1, C1, S1 = safeDiag!(Q1, Q2)
    # profiling
    # ans3 = @timed safeDiag!(Q1, Q2)
    # U1 = ans3[1][1]
    # V1 = ans3[1][2]
    # Z1 = ans3[1][3]
    # C1 = ans3[1][4]
    # S1 = ans3[1][5]

    # Step 4:
    # update U
    t = min(m, k+l)
    @views U[1:m,k+1:t] = U[1:m,k+1:t] * U1

    # Step 5:
    # update V
    @views V[1:p,1:l] = V[1:p,1:l] * V1

    # Step 6:
    # set W
    W = Z1' * R23

    # Step 7:
    # compute RQ decomposition of W
    row, col = size(W)
    W_ = copy(W)
    W, tau = LAPACK.gerqf!(W)
    Q3 = LAPACK.orgrq!(W, tau, length(tau))

    # Step 8:
    # update R13(a)
    @views A13 = A[1:k, n-l+1:n]*Q3'

    # Step 9:
    # update Q
    @views Q[1:n,n-l+1:n] = Q[1:n,n-l+1:n] * Q3'

    # form R
    R0 = zeros(T, k+l, n)
    @views R0[1:k, n-k-l+1:n-l+1] = A[1:k,n-k-l+1:n-l+1]
    @views R0[1:k, n-l+1:n] = A13
    @views R0[k+1:k+l, n-l+1:n]= W_ * Q3'

    # set C, S
    C = Matrix{T}(I, m, k+l)
    if m-k-l >= 0
        @views C[k+1:k+l,k+1:k+l] = C1
    else
        @views C[k+1:m,k+1:k+l] = C1
    end
    S = zeros(T, p, k+l)
    @views S[1:l, k+1:k+l] = S1
    GSVD(U, V, Q, C, S, k, l, R0)
    # profiling
    # return ans1[2], ans2[2], ans3[2]
end


function householderqr(R1, R2)
    r, c = size(R1)
    F = qr!([R1;R2])
    Q = Matrix(F.Q)
    Q1 = Q[1:r,:]
    Q2 = Q[r+1:end,:]
    return Q1, Q2, F.R
end

# end
