function Polynomials.derivative(p::TensorPolynomial{F,X,N,Y,M},z::Symbol) where {F,X,N,Y,M}
    z == X && return TensorPolynomial(derivative(p.px),p.py)
    z == Y && return TensorPolynomial(p.px,derivative(p.py))
    throw(ArgumentError("Z must be an indeterminate present in the field, but Z=$z was provided for a field with indeterminates X=$X and Y=$Y"))
end

# function Polynomials.derivative(p::BivariatePolynomial{F,X,N,Y,M},z::Symbol) where {F,X,N,Y,M}
#     z == X && return BivariatePolynomial(ImmutablePolynomial(derivative.(p.p.coeffs),:y))
#     z == Y && return BivariatePolynomial(derivative(p.p))
# end

function gradient(p::PolyScalarField{F,X,Y}) where {F,X,Y}
    dx = derivative(p,X)
    dy = derivative(p,Y)
    PolyVectorField(dx,dy)
end

function divergence(v::PolyVectorField)
    d1x = derivative(v.s1.px)
    d2y = derivative(v.s2.py)
    part1 = PolyScalarField(d1x,v.s1.py)
    part2 = PolyScalarField(v.s2.px,d2y)
    PolySum(part1,part2)
end

LinearAlgebra.dot(∇,v::PolyVectorField) = divergence(v)

laplacian(v::PolyScalarField) = divergence(gradient(v))

# This method should be removed in the next Polynomials update
Polynomials.derivative(p::ImmutablePolynomial{F,X,1}) where {F,X} = zero(p)



