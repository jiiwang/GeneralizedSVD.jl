# module Generalizedsvd

using LinearAlgebra
# using Main.CSD

import Base.show
include("preproc.jl")

struct GSVD{T,S} <: Factorization{T}
    U::S
    V::S
    Q::S
    D1::S
    D2::S
    k::Int
    l::Int
    R::S
    function GSVD{T,S}(U::AbstractMatrix{T}, V::AbstractMatrix{T}, Q::AbstractMatrix{T},
                                 D1::AbstractMatrix{T}, D2::AbstractMatrix{T}, k::Int, l::Int, R::AbstractMatrix{T}) where {T,S}
        new(U, V, Q, D1, D2, k, l, R)
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
        println(io, "\nR factor:")
        show(io, mime, F.R)
    end
end

function GSVD(U::AbstractMatrix{T}, V::AbstractMatrix{T}, Q::AbstractMatrix{T},
                        D1::AbstractMatrix{T}, D2::AbstractMatrix{T}, k::Int, l::Int, R::AbstractMatrix{T}) where T
    GSVD{T,typeof(U)}(U, V, Q, D1, D2, k, l, R)
end

# This function computes the generalized singular value
# decomposition (GSVD) of an m-by-n matrix A and p-by-n
#
#      A = U * D1 * R * Q',    B = V * D2 * R * Q'
#
#  where U, V and Q are orthogonal matrices. Let k+l be the effective numerical rank of the matrix (A',B')',
#  then R is a (k+l)-by-n matrix of structure [0 R0] where R0 is (k+l)-by-(k+l) and is nonsingular upper triangular matrix,
#  D1 and D2 are m-by-(k+l) and p-by-(k+l) "diagonal" matrices and
#  have one of the following structures:
#
#  If m >= k+l,
#
#         D1 =      k   l
#             k  (  I   0  )
#             l  (  0   C  )
#         m-k-l  (  0   0  )
#
#         D2 =      k   l
#             l  (  0   S  )
#           p-l  (  0   0  )
#
#
#  where
#
#       C = diag( alpha( k+1 ), ... , alpha( k+l ) ),
#       S = diag( beta( k+1 ), ..., beta( k+l ) ),
#       C**2 + S**2 = I
#
#
#  (2) If m < k+l,
#
#        D1 =         k   m-k  k+l-m
#                 k ( I    0     0  )
#               m-k ( 0    C     0  )
#
#         D2 = 	      k  m-k   k+l-m
#               m-k ( 0   S     0  )
#             k+l-m ( 0   0     I  )
#             p-k-l ( 0   0     0  )
#
#  where
#
#       C = diag( alpha(k+1), ... , alpha(m) ),
# 	    S = diag( beta(k+1), ..., beta(m) ),
#       C**2 + S**2 = I
#
#
# arguments:
# A: m by n matrix
# B: p by n matrix
# option: 0 or 1
#
# returns U, V, Z, alpha, beta, R, k, l if option = 0
# U, V, Z, D1, D2, R, k, l if option = 1

function gsvd(A, B)
# function gsvd(A, B, option)
    m, n = size(A)
    p = size(B)[1]

    # t_pre = 0.0
    # t_qr = 0.0
    # t_csd = 0.0

    # Step 1:
    # preprocess A, B
    # U, V, Q, k, l, A, B = dggsvp3!(A, B)
    # ans1 = @timed dggsvp3!(A, B)
    U, V, Q, k, l, A, B = preproc(A, B)
    # orthog_preproc(U,V,Q)

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
        R = [A[1:k,:];B[1:l,:]]
        if option == 0
            alpha = fill(0.0, k+l)
            for i = 1:m
                alpha[i] = 1.0
            end
            beta = fill(0.0, k+l)
            for i = m+1:k+l
                beta[i] = 1.0
            end
            return U, V, Q, alpha, beta, R, k, l
        else
            CS = [Matrix{Float64}(I, k+l, k+l);zeros(Float64, m+p-k-l, k+l)]
            C = @view CS[1:m, :]
            S = @view CS[m+1:m+p, :]
            # orthog_overall(U, V, Q)
            return U, V, Q, C, S, R, k, l
        end
    end

    # R1 = deepcopy(A23)
    # R2 = deepcopy(B13)
    # Q1, Q2, R23 = splitqr(A23, B13)

    Q1, Q2, R23 = householderqr(A23, B13)
    # ans2 = @timed householderqr(A23, B13)
    # t_qr = ans2[2]
    # Q1 = ans2[1][1]
    # Q2 = ans2[1][2]
    # R23 = ans2[1][3]
    # orthog_qr(Q1, Q2)

    # Alternatively, use Householder's method to do the qr of R1 and R2

    # Step 3:
    # CSD of Q1, Q2
    U1, V1, Z1, C1, S1 = csd(Q1, Q2, 1)
    # if option == 0
    #     alpha = fill(0.0, k+l)
    #     for i = 1:k
    #         alpha[i] = 1.0
    #     end
    #     beta = fill(0.0, k+l)
    #     U1, V1, Z1, alpha1, beta1 = csd(Q1, Q2, 0)
    #     if m-k-l >=0
    #         for i = k+1:l
    #             alpha[i] = alpha1[i]
    #             beta[i] = beta1[i]
    #         end
    #     else
    #         for i = 1:m-k
    #             alpha[k+i] = alpha1[i]
    #             beta[k+i] = beta1[i]
    #         end
    #         for i = m+1:k+l
    #             beta[i] = 1.0
    #         end
    #     end
    # else
    #     U1, V1, Z1, C1, S1 = csd(Q1, Q2, 1)
    #     # ans3 = @timed csd(Q1, Q2, 1)
    #     # t_csd = ans3[2]
    #     # U1 = ans3[1][1]
    #     # V1 = ans3[1][2]
    #     # Z1 = ans3[1][3]
    #     # C1 = ans3[1][4]
    #     # S1 = ans3[1][5]
    #     # orthog_csd(U1, V1, Z1)
    # end

    # Step 4:
    # update U
    t = min(m, k+l)
    @views U[1:m,k+1:t] = U[1:m,k+1:t] * U1

    # Step 5:
    # update V
    @views V[1:p,1:l] = V[1:p,1:l] * V1

    # Step 6:
    # set T
    T = Z1' * R23

    # Step 7:
    # compute RQ decomposition of T
    row, col = size(T)
    T_ = copy(T)
    T, tau = LAPACK.gerqf!(T)
    Q3 = LAPACK.orgrq!(T, tau, length(tau))

    # Step 8:
    # update R13(a)
    @views A13 = A[1:k, n-l+1:n]*Q3'

    # Step 9:
    # update Q
    @views Q[1:n,n-l+1:n] = Q[1:n,n-l+1:n] * Q3'

    # form R
    R = zeros(Float64, k+l, n)
    @views R[1:k, n-k-l+1:n-l+1] = A[1:k,n-k-l+1:n-l+1]
    @views R[1:k, n-l+1:n] = A13
    @views R[k+1:k+l, n-l+1:n]= T_ * Q3'

    # if option == 0
    #     return U, V, Q, alpha, beta, R, k, l
    # else
    #     # set C, S
    #     C = Matrix{Float64}(I, m, k+l)
    #     if m-k-l >= 0
    #         # for i in 1:l
    #         #     C[i+k,i+k] = C1[i,i]
    #         # end
    #         @views C[k+1:k+l,k+1:k+l] = C1
    #     else
    #         # for i in 1:m-k
    #         #     C[i+k,i+k] = C1[i,i]
    #         @views C[k+1:m,k+1:k+l] = C1
    #     end
    #     S = zeros(Float64, p, k+l)
    #     @views S[1:l, k+1:k+l] = S1
    #     # orthog_overall(U, V, Q)
    #     # return U, V, Q, C, S, R, k, l
    #     GSVD(U, V, Q, C, S, k, l, R)
    #     # return t_pre, t_qr, t_csd
    # end

    # set C, S
    C = Matrix{Float64}(I, m, k+l)
    if m-k-l >= 0
        # for i in 1:l
        #     C[i+k,i+k] = C1[i,i]
        # end
        @views C[k+1:k+l,k+1:k+l] = C1
    else
        # for i in 1:m-k
        #     C[i+k,i+k] = C1[i,i]
        @views C[k+1:m,k+1:k+l] = C1
    end
    S = zeros(Float64, p, k+l)
    @views S[1:l, k+1:k+l] = S1
    # orthog_overall(U, V, Q)
    # return U, V, Q, C, S, R, k, l
    GSVD(U, V, Q, C, S, k, l, R)
end


function householderqr(R1, R2)
    r, c = size(R1)
    F = qr([R1;R2])
    Q1 = Matrix(F.Q)[1:r,:]
    Q2 = Matrix(F.Q)[r+1:end,:]
    return Q1, Q2, F.R
end

# end
