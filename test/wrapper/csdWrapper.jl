const liblapack = Base.liblapack_name
import LinearAlgebra.BLAS.@blasfunc
import LinearAlgebra: BlasFloat, BlasInt, LAPACKException,
    DimensionMismatch, chkstride1, Givens

function dorcsd2by1!(X11::AbstractMatrix{Float64}, X21::AbstractMatrix{Float64})
    # chkstride1(X11, X21)
    jobu1 = 'Y'
    jobu2 = 'Y'
    jobv1t= 'Y'
    p, q = size(X11)
    m = p + size(X21, 1)
    ldx11 = max(1, p)
    ldx21 = max(1, m-p)
    theta = Vector{Float64}(undef, min(p, m-p, q, m-q))
    U1 = similar(X11, Float64, ldx11, p)
    U2 = similar(X21, Float64, ldx21, m-p)
    ldu1 = ldx11
    ldu2 = ldx21
    V1T = Matrix{Float64}(undef, q, q)
    ldv1t = max(1, q)

    work = Vector{Float64}(undef, 1)
    lwork = BlasInt(-1)
    iwork = Vector{BlasInt}(undef, m-min(p, m-p, q, m-q))
    info  = Ref{BlasInt}()


    for i in 1:2  # first call returns lwork as work[1]
        ccall((@blasfunc(dorcsd2by1_), liblapack), Cvoid,
        (Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
        Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
        Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
        Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
        Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}
        ),
        jobu1, jobu2, jobv1t,
        m, p, q,
        X11, ldx11, X21, ldx21,
        theta, U1, ldu1, U2, ldu2, V1T, ldv1t,
        work, lwork, iwork, info)
        # chklapackerror(info[])
        if i == 1
            lwork = BlasInt(real(work[1]))
            resize!(work, lwork)
        end
    end
    # return X11, X21, theta, U1, U2, V1T
    return theta, U1, U2, V1T
end
