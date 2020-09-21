function generator(row, col, rank)
    E = randn(rank, col)
    return B = randn(row ,rank)*E
end

function pre_gene(m, p, n, l, rA)
    B = generator(p, n, l)
    A = generator(m, n, rA)
    @time preproc(A, B)
end
