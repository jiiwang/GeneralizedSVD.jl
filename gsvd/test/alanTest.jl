Base.show(io::IO, x::Union{Float64,Float32}) = Base.Grisu._show(io, x, Base.Grisu.SHORTEST, 0, true, false)
m1,m2,n,r,rA,rB =2,2,3,2,1,2
v=[0,0,0]
trials = 200_000
for i = 1:trials

 E = randn(r,n)
 A = randn(m1,rA)*randn(rA,r)*E
 B = randn(m2,rB)*randn(rB,r)*E
 try
       U,V,Q,C,S,R =gsvd(A,B,1)
       v[size(C,2)] +=1
   catch
       println("Found an LAPACK Bug on trial #$i")
       show(A);show(B)
       println(v*100 ./trials)
   end

end
println("% of time, C is 2x2: ", v[2]*100 ./ trials)
println("% of time, C is 2x3: ", v[3]*100 ./ trials)
