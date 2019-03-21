using Pkg
Pkg.add("LinearAlgebra")
using LinearAlgebra

const liblapack = Base.liblapack_name
import LinearAlgebra.BLAS.@blasfunc
import LinearAlgebra: BlasFloat, BlasInt, LAPACKException,
    DimensionMismatch, chkstride1

function qsvd(A, B)
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
        A23 = A[k+1:k+l, n-l+1:n]
    # A23 is upper trapezoidal
    else
        A23 = A[k+1:m, n-l+1:n]
    end
    B13 = B[1:l, n-l+1:n]

    # R1 = copy(A23)
    # R2 = copy(B13)
    Q1, Q2, R23 = splitqr(A23, B13)

    # Step 3:
    # CSD of Q1, Q2
    U1, V1, Z1, C1, S1 = csd(Q1, Q2)

    # Step 4:
    # set C, S
    C = Matrix{Float64}(I, m, k+l)
    if m-k-l >= 0
        # for i in 1:l
        #     C[i+k,i+k] = C1[i,i]
        # end
        C[k+1:k+l,k+1:k+l] = C1
    else
        # for i in 1:m-k
        #     C[i+k,i+k] = C1[i,i]
        C[k+1:m,k+1:k+l] = C1
    end
    S = zeros(Float64, p, k+l)
    S[1:l, k+1:k+l] = S1

    # Step 5:
    # update U
    t = min(m, k+l)
    U[1:m,k+1:t] = U[1:m,k+1:t] * U1

    # Step 6:
    # update V
    V[1:p,1:l] = V[1:p,1:l] * V1

    # Step 7:
    # set T
    T = Z1' * R23

    # Step 8:
    # compute RQ decomposition of T
    row, col = size(T)
    T_ = copy(T)
    T, tau = LAPACK.gerqf!(T)
    Q3 = LAPACK.orgrq!(T, tau, length(tau))
    # R =  T_ * Q3'
    # R = UpperTriangular(T[1:row,col-row+1:col])
    # R * Q3 - T_

    # Step 9:
    # update R13(a)
    A13 = A[1:k, n-l+1:n]*Q3'

    # Step 10:
    # update Q
    Q[1:n,n-l+1:n] = Q[1:n,n-l+1:n] * Q3'

    # form R
    R = zeros(Float64, k+l, n)
    R[1:k, n-k-l+1:n-l+1] = A[1:k,n-k-l+1:n-l+1]
    R[1:k, n-l+1:n] = A13
    R[k+1:k+l, n-l+1:n]= T_ * Q3'

    return U, V, Q, C, S, R
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

# function preproc(A, B)
#     r = rank(B)
#     println(r)
#     F1 = qr(B, Val(true))
#     # F1.Q * F1.R * (F1.P)' - B
#     R = view(F1.R, 1:r,1:size(B,2))
#     # R = F1.R
#     m, n = size(R)
#     tau = zeros(Float64, min(m,n))
#     R, tau = LAPACK.gerqf!(R, tau)
#     # R, tau = LAPACK.gerqf!(R)
#     # R_ = copy(R)
#     # Q = LAPACK.orgrq!(R, tau, length(tau))
#     # F1.Q * B * F1.P * Q'
#     R
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
            if j <= m
                # eliminate element in R2 with elements in R1
                c, s, r = genGiv(R1[j, j], R2[i, j])
                R1[j, j] = r
                R2[i, j] = 0.0
                if j+1 <= n
                    R1[j, j+1:n], R2[i, j+1:n] = appGiv(R1[j, j+1:n], R2[i, j+1:n], c, s)
                end
            else
                # eliminate element in R2 with elements in R2
                c, s, r = genGiv(R2[j-m, j], R2[i, j])
                R2[j-m, j] = r
                R2[i, j] = 0.0
                if j+1 <= n
                    R2[j-m, j+1:n], R2[i, j+1:n] = appGiv(R2[j-m, j+1:n], R2[i, j+1:n], c, s)
                end
            end
            # update col in Q1 and Q2
            Q1[1:m, j], T[1:m] = appGiv(Q1[1:m, j], T[1:m], c, s)
            Q2[1:n, j], T[m+1:m+n] = appGiv(Q2[1:n, j], T[m+1:m+n], c, s)
        end
        if m+i <= n
           Q1[1:m,m+i] = T[1:m]
           Q2[1:n,m+i] = T[m+1:m+n]
        end
    end
    R = [R1; R2[1:n-m,:]]
    return Q1, Q2, R
end

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

function appGiv(v1, v2, c, s)
    u1 = c*v1 + s*v2
    u2 = c*v2 - s*v1
    return u1, u2
end
