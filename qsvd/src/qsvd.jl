using Pkg
Pkg.add("LinearAlgebra")
using LinearAlgebra

const liblapack = Base.liblapack_name
import LinearAlgebra.BLAS.@blasfunc
import LinearAlgebra: BlasFloat, BlasInt, LAPACKException,
    DimensionMismatch, chkstride1

function qsvd(A, B)
    # preprocess A & B
    dggsvp3!(A, B)
end

# function sggsvp3!(A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
#     ccall((@blasfunc($sggsvp3), liblapack), Cvoid, (Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
#     Int32, Int32, Int32,
#      ),)
#
#     JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA,
#     TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK,
#     WORK, LWORK, INFO
# end
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
            # resize!(iwork, iwork[1])
            # resize!(tau, tau[1])
        end
    end
    return U, V, Q, k, l
    # return U
end
#     end
# end

function preproc(A, B)
    r = rank(B)
    println(r)
    F1 = qr(B, Val(true))
    # F1.Q * F1.R * (F1.P)' - B
    R = view(F1.R, 1:r,1:size(B,2))
    # R = F1.R
    m, n = size(R)
    tau = zeros(Float64, min(m,n))
    R, tau = LAPACK.gerqf!(R, tau)
    # R, tau = LAPACK.gerqf!(R)
    # R_ = copy(R)
    # Q = LAPACK.orgrq!(R, tau, length(tau))
    # F1.Q * B * F1.P * Q'
    R
end
