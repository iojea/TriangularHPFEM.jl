module PolyFieldsTests

using Test
using TriangularhpFEM

p = TensorPolynomial((1.,2,5.),(2.,-1))
q = TensorPolynomial(((1.,2.),(-2,0.,1.)))
r = p+q
s = p*q
v = PolyVectorField([p,q])

x = HPPoint(1.,-1.)

w = v(x)

@test p(x) ≈ 24.
@test q(x) ≈ -3.
@test r(x) ≈ p(x)+q(x)
@test s(x) ≈ p(x)*q(x)
@test w[1] == p(x)
@test w[2] == q(x)

end; #module