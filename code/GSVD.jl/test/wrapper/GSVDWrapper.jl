using LinearAlgebra

const libblas = Base.libblas_name
const liblapack = Base.liblapack_name
import LinearAlgebra.BLAS.@blasfunc
import LinearAlgebra: BlasFloat, BlasInt, LAPACKException,
    DimensionMismatch, chkstride1, Givens

include("preprocWrapper.jl")

function dtgsja!(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64},
    U::AbstractMatrix{Float64}, V::AbstractMatrix{Float64}, Q::AbstractMatrix{Float64},
    k::Integer, l::Integer)
    # chkstride1(A, B)
    jobu = 'U'
    jobv = 'V'
    jobq = 'Q'
    m, n = size(A)
    p = size(B,1)
    tola = max(m, n) * opnorm(A, 1) * eps(Float64)
    tolb = max(p, n) * opnorm(B, 1) * eps(Float64)
    lda = max(1, stride(A, 2))
    ldb = max(1, stride(B, 2))
    ldu = max(1, m)
    ldv = max(1, p)
    ldq = max(1, n)
    alpha = Vector{Float64}(undef, n)
    beta = Vector{Float64}(undef, n)
    work = Vector{Float64}(undef, 2*n)
    ncycle  = Ref{BlasInt}()
    info  = Ref{BlasInt}()

    # for i in 1:2  # first call returns lwork as work[1]
    ccall((@blasfunc(dtgsja_), liblapack), Cvoid,
    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
    Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
    Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
    Ref{Float64}, Ref{Float64},
    Ptr{Float64}, Ptr{Float64},
    Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
    Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}
    ),
    jobu, jobv, jobq,
    m, p, n, k, l,
    A, lda, B, ldb,
    tola, tolb,
    alpha, beta,
    U, ldu, V, ldv, Q, ldq,
    work, ncycle, info)
    # end
    return A, B, alpha, beta, U, V, Q, work, ncycle, info
end

function gsvdwrapper(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64})
    m, n = size(A)
    p = size(B)[1]
    # U, V, Q, k, l, A, B = dggsvp3!(A, B)
    # profiling
    ans1 = @timed dggsvp3!(A, B)
    U = ans1[1][1]
    V = ans1[1][2]
    Q = ans1[1][3]
    k = ans1[1][4]
    l = ans1[1][5]
    A = ans1[1][6]
    B = ans1[1][7]

    # A, B, alpha, beta, U, V, Q, work, ncycle, info = dtgsja!(A, B, U, V, Q, k, l)
    # profiling
    ans2 = @timed dtgsja!(A, B, U, V, Q, k, l)
    A = ans2[1][1]
    B = ans2[1][2]
    alpha = ans2[1][3]
    beta = ans2[1][4]
    U = ans2[1][5]
    V = ans2[1][6]
    Q = ans2[1][7]

    R0 = zeros(k+l, n)
    D1 = Matrix(1.0I, m, k+l)
    D2 = zeros(p, k+l)
    if m-k-l >= 0
        @views R0[:, n-k-l+1:n]  = A[1:k+l, n-k-l+1:n]
        @views D1[k+1:k+l, k+1:k+l] = Diagonal(alpha[k+1:k+l])
        @views D2[1:l, k+1:k+l] = Diagonal(beta[k+1:k+l])
    else
        # println("$k $l")
        @views R0[1:m, n-k-l+1:n] = A[1:m, n-k-l+1:n]
        @views R0[m+1:k+l, m+n-k-l+1:n] = B[m-k+1:l, n+m-k-l+1:n]
        @views D1[k+1:m, k+1:m] = Diagonal(alpha[k+1:m])
        @views D2[1:m-k, k+1:m] = Diagonal(beta[k+1:m])
        @views D2[m-k+1:l, m+1:k+l] = Matrix(1.0I, k+l-m, k+l-m)
    end
    # return U, V, Q, D1, D2, k, l, R0
    # profiling
    return ans1[2], ans2[2]
end

# function compute_dot(DX::Vector{Float64}, DY::Vector{Float64})
#     @assert length(DX) == length(DY)
#     n = length(DX)
#     incx = incy = 1
#     ccall((@blasfunc(ddot_), libblas), Float64,
#       (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
#       Ptr{Float64}, Ref{BlasInt}),
#       n, DX, incx, DY, incy)
# end
#
# function compute_dot_julia(dx::Vector{Float64}, dy::Vector{Float64})
#     # ddot = 0.0
#     dtemp = 0.0
#     n = length(dx)
#     m = n % 5
#     if m != 0
#         for i = 1:m
#             dtemp = dtemp + dx[i]*dy[i]
#         end
#         if n < 5
#             # ddot = dtemp
#             return dtemp
#         end
#     end
#     mp1 = m + 1
#     for i = mp1:5:n-5
#         dtemp = dtemp + dx[i]*dy[i] + dx[i+1]*dy[i+1] +
#               dx[i+2]*dy[i+2] + dx[i+3]*dy[i+3] + dx[i+4]*dy[i+4]
#     end
#     # ddot = dtemp
#     # return ddot
#     return dtemp
# end
