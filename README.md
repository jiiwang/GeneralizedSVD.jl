# GeneralizedSVD

*A Julia program to compute the generalized singular value decomposition (GSVD).*

## Installation
To install `GeneralizedSVD`,
from the Julia REPL, type `]` to enter the Pkg REPL mode and run:
```
pkg> add GeneralizedSVD
```

or using the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("GeneralizedSVD")
```

## Example
```julia
julia> using GeneralizedSVD

julia> A = [1. 2 3 0; 5 4 2 1; 0 3 5 2; 2 1 3 3; 2 0 5 3];
julia> B = [1. 0 3 -1; -2 5 0 1; 4 2 -1 2];

julia> F = gsvd(A, B);

julia> F.C
5×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 4 stored entries:
  [1, 1]  =  1.0
  [2, 2]  =  0.894685
  [3, 3]  =  0.600408
  [4, 4]  =  0.27751

julia> F.S
3×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 3 stored entries:
  [1, 2]  =  0.446698
  [2, 3]  =  0.799694
  [3, 4]  =  0.960723

julia> F.R
4×4 Array{Float64,2}:
 5.74065  -7.07986   0.125979  -0.316232
 0.0      -7.96103  -2.11852   -2.98601
 0.0       0.0       5.72211   -0.43623
 0.0       0.0       0.0        5.66474

julia> [A; B] ≈ [F.U*F.C; F.V*F.S]*F.R*F.Q'
true
```

## Documentation

- [docs-url] &mdash; **documentation of the most recently tagged version.**

## Project Status

The package is tested against, and being developed for, Julia `1.3` on macOS and supports `Float` type.

## Contact
+ Ji Wang: jiiwang@ucdavis.edu
+ Zhaojun Bai: zbai@ucdavis.edu

[docs-url]:https://jiiwang.github.io/GeneralizedSVD.jl/
