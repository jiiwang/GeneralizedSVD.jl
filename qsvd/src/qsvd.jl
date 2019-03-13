using Pkg
Pkg.add("LinearAlgebra")
using LinearAlgebra


function qsvd(A, B)
    # preprocess A & B
    sggsvp3!(A, B)
    # ccall((:clock, "msvcrt"), Int32, ())
end

function sggsvp3!(A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
    ccall((@blasfunc($sggsvp3), liblapack), Cvoid, (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, 
    Int32, Int32, Int32,
     ),)

    JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA,
    TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK,
    WORK, LWORK, INFO
end
