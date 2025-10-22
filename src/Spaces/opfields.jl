struct PolySum{F,X,Y,F1<:PolyField{F,X,Y},F2<:PolyField{F,X,Y}} <: PolyScalarField{F,X,Y}
    left::F1
    right::F2
end

(p::PolySum)(x,y) = p.left(x,y)+p.right(x,y)
(p::PolySum)(x) = p.left(x) + p.right(x)

Base.:+(p::P,q::Q) where {F,X,Y,P<:PolyField{F,X,Y},Q<:PolyField{F,X,Y}} = PolySum(p,q)

LinearAlgebra.dot(p::PolyVectorField{F,X,Y},q::PolyVectorField{F,X,Y})  where {F,X,Y} = p.s1*q.s1 + p.s2*q.s2

Base.:*(ps::PolySum,p::TensorPolynomial) = ps.left*p + ps.right*p
Base.:*(p::TensorPolynomial,ps::PolySum) = ps*p
Base.:*(ps::PolySum,p::PolyVectorField) = PolyVectorField(ps*p.s1,ps*p.s2)
Base.:*(p::PolyVectorField,ps::PolySum) = ps*p
Base.:*(ps::PolySum,qs::PolySum) = ps.left*qs.left + ps.right*qs.left + ps.right*qs.left + ps.right*qs.right







