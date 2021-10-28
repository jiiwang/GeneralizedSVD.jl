module GSVD

import Base: show,getproperty
using LinearAlgebra
using Base: promote_typeof

export GenSVD, gsvd, gsvd!

include("preproc.jl")
include("householderqr.jl")
include("csd.jl")

# Immutable composite type GenSVD
struct GenSVD{T,S} <: Factorization{T}
    U::S
    V::S
    Q::S
    alpha::Vector
    beta::Vector
    k::Int
    l::Int
    R1::S

    # Inner constructor
    function GenSVD{T,S}(U::AbstractMatrix{T}, V::AbstractMatrix{T},
                       Q::AbstractMatrix{T}, alpha::Vector,
                       beta::Vector, k::Int, l::Int,
                       R1::AbstractMatrix{T}) where {T,S}
        new(U, V, Q, alpha, beta, k, l, R1)
    end
end

# Outer constructor
function GenSVD(U::AbstractMatrix{T}, V::AbstractMatrix{T},
              Q::AbstractMatrix{T}, alpha::Vector,
              beta::Vector, k::Int, l::Int,
              R1::AbstractMatrix{T}) where T
    GenSVD{T,typeof(U)}(U, V, Q, alpha, beta, k, l, R1)
end

# iteration for destructuring into components
Base.iterate(S::GenSVD) = (S.U, Val(:V))
Base.iterate(S::GenSVD, ::Val{:V}) = (S.V, Val(:Q))
Base.iterate(S::GenSVD, ::Val{:Q}) = (S.Q, Val(:C))
Base.iterate(S::GenSVD, ::Val{:C}) = (S.C, Val(:S))
Base.iterate(S::GenSVD, ::Val{:S}) = (S.S, Val(:R))
Base.iterate(S::GenSVD, ::Val{:R}) = (S.R, Val(:done))
Base.iterate(S::GenSVD, ::Val{:done}) = nothing

# Explicitely construct C, S, vals and R
@inline function getproperty(F::GenSVD{T}, d::Symbol) where T
    Fa = getfield(F, :alpha)
    Fb = getfield(F, :beta)
    Fk = getfield(F, :k)
    Fl = getfield(F, :l)
    FU = getfield(F, :U)
    FV = getfield(F, :V)
    FQ = getfield(F, :Q)
    FR = getfield(F, :R1)
    if d === :vals
        return Fa[1:Fk + Fl] ./ Fb[1:Fk + Fl]
    elseif d === :C
        m = size(FU, 1)
        if m - Fk - Fl >= 0
            return [Matrix{T}(I, Fk, Fk)  zeros(T, Fk, Fl)            ;
                    zeros(T, Fl, Fk)      Diagonal(Fa[Fk + 1:Fk + Fl]);
                    zeros(T, m - Fk - Fl, Fk + Fl)                    ]
        else
            return [Matrix{T}(I, m, Fk) [zeros(T, Fk, m - Fk); Diagonal(Fa[Fk + 1:m])] zeros(T, m, Fk + Fl - m)]
        end
    elseif d === :S
        m = size(FU, 1)
        p = size(FV, 1)
        if m - Fk - Fl >= 0
            return [zeros(T, Fl, Fk) Diagonal(Fb[Fk + 1:Fk + Fl]); zeros(T, p - Fl, Fk + Fl)]
        else
            return [zeros(T, p, Fk) [Diagonal(Fb[Fk + 1:m]); zeros(T, Fk + p - m, m - Fk)] [zeros(T, m - Fk, Fk + Fl - m); Matrix{T}(I, Fk + p - m, Fk + Fl - m)]]
        end
    elseif d === :R
        n = size(FQ, 1)
        return [zeros(T, Fk + Fl, n - Fk - Fl) FR]
    else
        getfield(F, d)
    end
end

Base.propertynames(F::GenSVD) =
    (:vals, :C, :S, :R, fieldnames(typeof(F))...)

# show() is used for neat display
function show(io::IO, mime::MIME{Symbol("text/plain")},
              F::GenSVD{<:Any,<:AbstractArray})
    summary(io, F); println(io)
    println(io, "U factor:")
    show(io, mime, F.U)
    println(io, "\n V factor:")
    show(io, mime, F.V)
    println(io, "\n Q factor:")
    show(io, mime, F.Q)
    println(io, "\n C factor:")
    show(io, mime, F.C)
    println(io, "\n S factor:")
    show(io, mime, F.S)
    println(io, "\n R factor:")
    show(io, mime, F.R)
    println(io, "\n Generalized singular values:")
    show(io, mime, F.vals)
    println(io, "\n k factor:")
    show(io, mime, F.k)
    println(io, "\n l factor:")
    show(io, mime, F.l)
    println(io, "\n numerical rank of the matrix [A; B]:")
    show(io, mime, F.k + F.l)
end

# promotion type to use for eigenvalues of a Matrix{T}
eigtype(T) = promote_type(Float32, typeof(zero(T)/sqrt(abs2(one(T)))))

copy_oftype(A::AbstractArray{T}, ::Type{T}) where {T} = copy(A)
copy_oftype(A::AbstractArray{T,N}, ::Type{S}) where {T,N,S} = convert(AbstractArray{S,N}, A)

"""
    gsvd(A, B) -> GenSVD

Compute the generalized SVD of an m-by-n matrix `A` and a p-by-n matrix `B`,
returning a `GenSVD` factorization object `F` such that

```math
     A = F.U * F.C * F.R * F.Q',    B = F.V * F.S * F.R * F.Q'
```

# Examples
```jldoctest
julia> A = randn(3, 2); B = randn(4, 2);

julia> F = gsvd(A, B);

julia> [A; B] ≈ [F.U*F.C; F.V*F.S]*F.R*F.Q'
true
```
"""
function gsvd(A::StridedMatrix{T}, B::StridedMatrix{T}) where T<:AbstractFloat
    return gsvd!(copy(A),copy(B))
end

function gsvd(A::StridedMatrix{TA}, B::StridedMatrix{TB}) where {TA,TB}
    S = promote_type(eigtype(TA),TB)
    return gsvd!(copy_oftype(A, S), copy_oftype(B, S))
end

"""
    gsvd!(A, B) -> GenSVD
`gsvd!` is the same as [`gsvd`](@ref), but modifies the arguments
`A` and `B` in-place, instead of making copies.
"""
function gsvd!(A::StridedMatrix{T}, B::StridedMatrix{T}) where T<:AbstractFloat
    m, n = size(A)
    p = size(B)[1]

    if n != size(B)[2]
        throw(ArgumentError("A's column number must equal that of B"))
    end

    # Step 1:
    # preprocess A, B
    U, V, Q, k, l, A, B = preproc!(A, B)

    # Early exit if A23 doesn't exist
    if l == 0
        return GenSVD(U, V, Q, ones(T, k), zeros(T, k), k, 0, A[1:k+l,n-k-l+1:n])
    end

    if m == k
        return GenSVD(U, V, Q, [ones(T, k); zeros(T, l)], [zeros(T, k); ones(T, l)], k, l, [A[1:k,n-k-l+1:n];B[1:l,n-k-l+1:n]])
    end

    # A23 is upper triangular
    if m-k-l >= 0
        @views A23 = A[k+1:k+l, n-l+1:n]
    # A23 is upper trapezoidal
    else
        @views A23 = A[k+1:m, n-l+1:n]
    end
    @views B13 = B[1:l, n-l+1:n]

    # Step 2:
    # QR decomposition of A23 and B13
    Q1, Q2, R23 = householderqr(A23, B13)

    # Step 3:
    # CSD of Q1, Q2
    U1, V1, Z1, cosine, sine = csd2by1!(Q1, Q2)

    # Step 4:
    # set alpha, beta
    alpha = ones(T, k+l)
    @views alpha[k+1:k+l] = cosine
    beta = zeros(T, k+l)
    @views beta[k+1:k+l] = sine

    # Step 5:
    # update U
    t = min(m, k+l)
    U[1:m,k+1:t] = @views U[1:m,k+1:t] * U1

    # Step 6:
    # update V
    V[1:p,1:l] = @views V[1:p,1:l] * V1

    # Step 7:
    # set W
    W = Z1' * R23

    # Step 8:
    # compute RQ decomposition of W
    W_ = copy(W)
    W, tau = LAPACK.gerqf!(W)
    Q3 = LAPACK.orgrq!(W, tau, length(tau))

    # Step 9:
    # update R13(a)
    @views A13 = A[1:k, n-l+1:n]*Q3'
    # form R1
    R1 = zeros(T, k+l, k+l)
    @views R1[1:k, 1:k+1] = A[1:k,n-k-l+1:n-l+1]
    @views R1[1:k, k+1:k+l] = A13
    @views R1[k+1:k+l, k+1:k+l]= W_ * Q3'

    # Step 10:
    # update Q
    @views Q[1:n,n-l+1:n] = Q[1:n,n-l+1:n] * Q3'

    GenSVD(U, V, Q, alpha, beta, k, l, triu!(R1))
end

end
