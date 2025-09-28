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


struct OperationField{F<:Function,F1<:Function,F2<:Function} <:Function
    op::F
    args::Tuple{F1,F2}
end
(of::OperationField)(x) = of.F(of.args1(x),of.args2(x))

function (v::PolyVectorField)(x::T) where T<:AbstractVector
    T(v.s1(x),v.s2(x))
end

struct OperationField{F<:Function,F1<:Function,F2<:Function} <:Function
    op::F
    args::Tuple{F1,F2}
end
(of::OperationField)(x) = of.F(of.args1(x),of.args2(x))




function Polynomials.:derivative(p::PolyScalarField{F,X,N,Y,M},z::Symbol) where {F,X,N,Y,M}
    z == X && return PolyScalarField(derivative(p.px),p.py)
    z == Y && return PolyScalarField(p.px,derivative(p.py))
    throw(ArgumentError("Z must be an indeterminate present in the field, but Z=$z was provided for a field with indeterminates X=$X and Y=$Y"))
end

function gradient(p::PolyScalarField{F,X,N,Y,M}) where {F,X,N,Y,M}
    dx = derivative(p,X)
    dy = derivative(p,Y)
    PolyVectorField(dx,dy)
end

function div(v::PolyVectorField)
    d1x = derivative(v.s1.px)
    d2y = derivative(v.s2.py)
    part1 = PolyScalarField(d1x,v.s1.py)
    part2 = PolyScalarField(v.s2.px,d2y)
    OperationField(+,(part1,part2))
end

