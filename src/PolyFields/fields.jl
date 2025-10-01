abstract type PolyField <:Function end


###########################
#     SCALAR FIELDS
###########################



indeterminate(::AbstractPolynomial{T,X}) where {T,X} = X
"""
```
    julia> ϕ = PolyScalarField((2,1.,3.),(3.,2.))
    julia> ϕ(1,-1)
        6.0
    julia> ϕ([1,-1])
        6.0
```
Defines a scalar field `ϕ:ℝ²→ℝ` given by a product of a polynomial in `x` and  a polynomial in `y`. In this case: `(2.0 + 1.0x + 3.0x²)(3.0+2.0y)`.

The indeterminates can be specified:
```
    julia> ϕ= PolyScalarField((2,1.,3.),(3.,2.),:z,:ξ)
```
produces the same polynomials, but on the variables `(z,ξ)`. 
"""
struct PolyScalarField{F,X,N,Y,M} <: PolyField
    px::ImmutablePolynomial{F,X,N}
    py::ImmutablePolynomial{F,Y,M}
end

function PolyScalarField(t1::Tuple,t2::Tuple,X,Y)
    t1 = promote(t1...)
    t2 = promote(t2...)
    F = promote_type(eltype(t1),eltype(t2))
    N = length(t1)
    M = length(t2)
    p1 = ImmutablePolynomial(NTuple{N,F}(convert.(F,t1)),X)
    p2 = ImmutablePolynomial(NTuple{M,F}(convert.(F,t2)),Y)
    PolyScalarField{F,X,N,Y,M}(p1,p2)
end
PolyScalarField(t1::Tuple,t2::Tuple) =  PolyScalarField(t1,t2,:x,:y)
(s::PolyScalarField)(x,y) = s.px(x)*s.py(y)
(s::PolyScalarField)(x::T) where T<:AbstractVector = s.px(x[1])*s.py(x[2])


indeterminates(p::PolyScalarField) = (indeterminate(p.px),indeterminate(p.py))

Base.:*(a::Number,p::PolyScalarField) = PolyScalarField(a*p.px,p.py)

function Base.:*(p::AbstractPolynomial,q::PolyScalarField)
    if indeterminate(p)===indeterminate(q.px)
        PolyScalarField(p*q.px,q.py)
    elseif indeterminate(p)===indeterminate(q.py)
        PolyScalarField(q.px,p*q.py)
    else
        throw(ArgumentError("Indeterminates does not match. Polynomial p has indeterminate $(indeterminate(p)), whereas PolyScalarField q has indeterminates $(indeterminate(q.px)),$(indeterminate(q.py))."))
    end
end
Base.:*(q::PolyScalarField,p::AbstractPolynomial) = p*q

function Base.:*(p::P,q::Q) where {P<:PolyScalarField,Q<:PolyScalarField}
    PolyScalarField(p.px*q.px,p.py*q.py)
end


###############################
#       VECTOR FIELDS
###############################
struct PolyVectorField{F,X,N1,N2,Y,M1,M2} <: PolyField
    s1::PolyScalarField{F,X,N1,Y,M1}
    s2::PolyScalarField{F,X,N2,Y,M2}
end

(v::PolyVectorField)(::Type{T},x,y) where T= T(v.s1(x,y),v.s2(x,y))
(v::PolyVectorField)(x,y) = [v.s1(x,y),v.s2(x,y)]
(v::PolyVectorField)(x::T) where T<:AbstractVector = v(T,x[1],x[2])

indeterminates(v::PolyVectorField) = indeterminates(v.s1)

function Base.:*(p::PolyScalarField,v::PolyVectorField)
    issubset(indeterminates(p),indeterminates(v)) || throw(ArgumentError("Fields have different indeterminates"))
    PolyVectorField(p*v.s1,p*v.s2)
end
function Base.:*(p::AbstractPolynomial,v::PolyVectorField)
    indeterminate(p) in indeterminates(v) || throw(ArgumentError("Indeterminates does not match."))
    PolyVectorField(p*v.s1,p*v.s2)    
end
Base.:*(a::Number,v::PolyVectorField) = PolyVectorField(a*v.s1,a*v.s2)

