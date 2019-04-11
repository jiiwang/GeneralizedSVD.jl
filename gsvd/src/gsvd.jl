using Pkg
Pkg.add("LinearAlgebra")
using LinearAlgebra

const liblapack = Base.liblapack_name
import LinearAlgebra.BLAS.@blasfunc
import LinearAlgebra: BlasFloat, BlasInt, LAPACKException,
    DimensionMismatch, chkstride1, Givens

# This function computes the quotient (or generalized) singular value
# decomposition (QSVD/GSVD) of an m-by-n matrix A and p-by-n
# matrix B:
#
#      A = U * D1 * R * Q',    B = V * D2 * R * Q'
#
#  where U, V and Q are orthogonal matrices. Let k+l be the effective numerical rank of the matrix (A',B')',
#  then R is a (K+L)-by-n nonsingular upper triangular matrix,
#  D1 and D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and
#  have one of the following structures:
#
#  If M >= K+L,
#
#         D1 =      K   L
#             K  (  I   0  )
#             L  (  0   C  )
#         M-K-L  (  0   0  )
#
#         D2 =      K   L
#             L  (  0   S  )
#           P-L  (  0   0  )
#
#
#  where
#
#       C = diag( ALPHA( K+1 ), ... , ALPHA( K+L ) ),
#       S = diag( BETA( K+1 ), ..., BETA( K+L ) ),
#       C**2 + S**2 = I
#
#
#  (2) If M < K+L,
#
#        D1 =         K   M-K  K+L-M
#                 K ( I    0     0  )
#               M-K ( 0    C     0  )
#
#         D2 = 	      K  M-K   K+L-M
#               M-K ( 0   S     0  )
#             K+L-M ( 0   0     I  )
#             P-K-L ( 0   0     0  )
#
#  where
#
#       C = diag( ALPHA(K+1), ... , ALPHA(M) ),
# 	    S = diag( BETA(K+1), ..., BETA(M) ),
#       C**2 + S**2 = I
#
# argument
# A: m by n matrix
# B: p by n matrix
#
# returns U, V, Z, alpha, beta, R, k, l
function gsvd(A, B, option)
    m, n = size(A)
    p = size(B)[1]

    # Step 1:
    # preprocess A, B
    U, V, Q, k, l, A, B = dggsvp3!(A, B)
    # rank([A' B']') == k+l

    # Step 2:
    # QR in a split fashion
    # A23 is upper triangular
    if m-k-l >= 0
        @views A23 = A[k+1:k+l, n-l+1:n]
    # A23 is upper trapezoidal
    else
        @views A23 = A[k+1:m, n-l+1:n]
    end
    @views B13 = B[1:l, n-l+1:n]

    # R1 = copy(A23)
    # R2 = copy(B13)
    Q1, Q2, R23 = splitqr(A23, B13)

    # Step 3:
    # CSD of Q1, Q2
    if option == 0
        alpha = fill(0.0, n)
        for i = 1:k
            alpha[i] = 1.0
        end
        beta = fill(0.0, n)
        U1, V1, Z1, alpha1, beta1 = csd(Q1, Q2, 0)
        if m-k-l >=0
            for i = 1:l
                alpha[k+i] = alpha1[i]
                beta[k+i] = beta1[i]
            end
        else
            for i = 1:m-k
                alpha[k+i] = alpha1[i]
                beta[k+i] = beta1[i]
            end
            for i = m+1:k+l
                beta[i] = 1.0
            end
        end
    else
        U1, V1, Z1, C1, S1 = csd(Q1, Q2, 1)
    end

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

    if option == 0
        return U, V, Q, alpha, beta, R, k, l
    else
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
        return U, V, Q, C, S, R, k, l
    end
end


# for (sggsvp3, elty) in (:sggsvp3_,:Float64)
    # @eval begin
        # function sggsvp3!(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64}, tau::AbstractVector{Float64})
function dggsvp3!(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64})
    chkstride1(A, B)
    jobu = 'U'
    jobv = 'V'
    jobq = 'Q'
    m, n = size(A)
    p = size(B,1)
    # if BlasInt(size(B, 2)) != n
    #     throw(DimensionMismatch("B has second dimension $(size(B,2)) but needs $n"))
    # end
    tola = max(m, n) * norm(A, 1) * eps(Float64)
    tolb = max(p, n) * norm(B, 1) * eps(Float64)
    k = Ref{BlasInt}()
    l = Ref{BlasInt}()
    ldu = max(1, m)
    U = similar(A, Float64, ldu, m)
    ldv = max(1, p)
    V = similar(A, Float64, ldv, p)
    ldq = max(1, n)
    Q = similar(A, Float64, ldq, n)
    iwork = Vector{BlasInt}(undef, n)
    tau = Vector{Float64}(undef, n)
    work = Vector{Float64}(undef, 1)
    lwork = BlasInt(-1)
    info  = Ref{BlasInt}()

    for i in 1:2  # first call returns lwork as work[1]
        ccall((@blasfunc(dggsvp3_), liblapack), Cvoid,
        (Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
        Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
        Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt},
        Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
        Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}
        ),
        jobu, jobv, jobq,
        m, p, n, A, max(1, stride(A, 2)), B, max(1, stride(B, 2)),
        tola, tolb, k, l,
        U, ldu, V, ldv, Q, ldq,
        iwork, tau, work, lwork, info)
        # chklapackerror(info[])
        if i == 1
            lwork = BlasInt(real(work[1]))
            resize!(work, lwork)
        end
    end
    return U, V, Q, k[], l[], A, B
end
#     end
# end

function splitqr(R1, R2)
    m, n = size(R1)

    # Initialize Q1, Q2 so that [Q1' Q2']' = I'
    Q1 = Matrix{Float64}(I, m, n)
    Q2 = zeros(Float64, n, n)
    for i = 1:n-m
        Q2[i, i+m] = 1.0
    end

    for i = 1:n
        # Initialize T: m+n vector
        T = fill(0.0, m+n)
        T[m+i] = 1.0
        k = min(n, m-1+i)
        for j = 1:k
            flag = 0
            if j <= m
                # eliminate element in R2 with elements in R1
                c, s, r = genGiv(R1[j, j], R2[i, j])
                # if j != i
                #     flag = 1
                    # c, s, r = givensAlgorithm(R1[j, j], R2[i, j])
                    R1[j, j] = r
                    R2[i, j] = 0.0
                    if j+1 <= n
                        # @views Ra = R1[j, j+1:n]
                        # @views Rb = R2[i, j+1:n]
                        R1[j, j+1:n], R2[i, j+1:n] = appGiv(R1[j, j+1:n], R2[i, j+1:n], c, s)
                        # R1[j, j+1:n] = G * R1[j, j+1:n]
                        # R2[i, j+1:n] = G * R2[i, j+1:n]
                    end
                # end
            else
                # eliminate element in R2 with elements in R2
                c, s, r = genGiv(R2[j-m, j], R2[i, j])
                # if j-m != i
                    # flag = 1
                    # c, s, r = givensAlgorithm(R2[j-m, j], R2[i, j])
                    R2[j-m, j] = r
                    R2[i, j] = 0.0
                    if j+1 <= n
                        # @views Ra = R2[j-m, j+1:n]
                        # @views Rb = R2[i, j+1:n]
                        R2[j-m, j+1:n], R2[i, j+1:n] = appGiv(R2[j-m, j+1:n], R2[i, j+1:n], c, s)
                        # @views R2[j-m, j+1:n] = [c s;-s c]*[R2[j-m, j+1:n], R2[i, j+1:n]][1]
                        # @views R2[i, j+1:n] = [c s;-s c]*[R2[j-m, j+1:n], R2[i, j+1:n]][2]
                        # R2[j-m, j+1:n] = G * R2[j-m, j+1:n]
                        # R2[i, j+1:n] = G * R2[i, j+1:n]
                    end
                # end
            end
            # update col in Q1 and Q2
            # @views Ra = R1[j, j+1:n]
            # @views Rb = R2[i, j+1:n]
            Q1[1:m, j], T[1:m] = appGiv(Q1[1:m, j], T[1:m], c, s)
            Q2[1:n, j], T[m+1:m+n] = appGiv(Q2[1:n, j], T[m+1:m+n], c, s)
            # if flag == 1
                # Q1[1:m, j] = G * Q1[1:m, j]
                # T[1:m] = G * T[1:m]
                # Q2[1:n, j] = G * Q2[1:n, j]
                # T[m+1:m+n] = G * T[m+1:m+n]
                # flag = 0
            # end
        end
        if m+i <= n
           @views Q1[1:m,m+i] = T[1:m]
           @views Q2[1:n,m+i] = T[m+1:m+n]
        end
    end
    @views R = [R1; R2[1:n-m,:]]
    return Q1, Q2, R
end

# This function generates a 2 by 2 Givens rotation matrix
function genGiv(a, b)
    c = 0.0
    s = 1.0
    r = b
    if a != 0
        if abs(a) > abs(b)
            t = b/a
            u = sqrt(1 + t^2)
            c = 1/u
            s = c*t
            r = a*u
        else
            t = a/b
            u = sqrt(1 + t^2)
            s = 1/u
            c = s*t
            r = b*u
        end
    end
    return c, s, r
end

# function appGiv(v1, v2, c, s)
#     u1 = c*v1 + s*v2
#     u2 = c*v2 - s*v1
#     return u1, u2
# end

# This function applies a 2 by 2 Givens rotation matrix
function appGiv(v1, v2, c, s)
    return c*v1 + s*v2, c*v2 - s*v1
end