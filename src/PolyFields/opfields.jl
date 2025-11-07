struct PolySum{F,X,Y} <: PolyScalarField{F,X,Y}
    left::PolyScalarField{F,X,Y}
    right::PolyScalarField{F,X,Y}
end

(p::PolySum)(x,y) = p.left(x,y)+p.right(x,y)
(p::PolySum)(x) = p.left(x) + p.right(x)

function Base.:+(p::P,q::Q) where {F,X,Y,P<:PolyField{F,X,Y},Q<:PolyField{F,X,Y}}
    if p==zero(p)
        q
    elseif q == zero(q)
        p
    else
        PolySum(p,q)
    end
end

function Base.zero(::PolySum{F,X,Y}) where {F,X,Y}
    zero(TensorPolynomial{F,X,Y})
end
#LinearAlgebra.dot(p::PolyVectorField{F,X,Y},q::PolyVectorField{F,X,Y})  where {F,X,Y} = p.s1*q.s1 + p.s2*q.s2

Base.:*(ps::PolySum,p::TensorPolynomial) = ps.left*p + ps.right*p
Base.:*(p::TensorPolynomial,ps::PolySum) = ps*p
Base.:*(ps::PolySum,p::PolyVectorField) = PolyVectorField(ps*p.s1,ps*p.s2)
Base.:*(p::PolyVectorField,ps::PolySum) = ps*p
Base.:*(ps::PolySum,qs::PolySum) = ps.left*qs.left + ps.right*qs.left + ps.right*qs.left + ps.right*qs.right

