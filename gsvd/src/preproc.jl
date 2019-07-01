using LinearAlgebra

function preproc(A, B)
    # Step 1: compute QR decomposition
    # with col pivoting of matrix B
    F = qr(B, Val(true))
end
