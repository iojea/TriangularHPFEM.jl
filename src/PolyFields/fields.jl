abstract type PolyField{F,X,Y} <:Function end

###########################################################################
#################              SCALAR FIELDS              #################
###########################################################################
#
abstract type PolyScalarField{F,X,N,Y,M} <: PolyField{F,X,Y} end
indeterminate(::AbstractPolynomial{T,X}) where {T,X} = X
indeterminates(::PolyField{F,X,Y}) where {F,X,Y} = (X,Y)
# 
######################################
####      BivariatePolynomial     ####
######################################

_zerofill(t::NTuple{K,T},N) where {K,T} = K>=N ? t : (t...,zeros(T,N-K)...)

"""
```
    BivariatePolynomial(t::Tuple)
```
Creates a bivariate polynomial with the coefficients stored in  `t`. `t` should be a tuple of tuples. Each sub-tuple defines a polynomial on the variable `x`. This polynomials are the coefficients of a polynomial on `y`, thus forming a bivariate polynomial.

""" 
struct BivariatePolynomial{F,X,N,Y,M} <: PolyScalarField{F,X,N,Y,M}
    p::ImmutablePolynomial{ImmutablePolynomial{F,X,N},Y,M}
end
function BivariatePolynomial{X,Y}(t::Tuple) where {X,Y}
        all(typeof(u)<:Tuple for u in t) || throw(ArgumentError("A tuple of tuples is expected"))
        N = maximum(length.(t))
        M = length(t)
        t = Tuple(promote(u...) for u in t)
        F = promote_type(eltype.(t)...)
        q = Tuple(ImmutablePolynomial(_zerofill(F.(u),N),:x) for u in t)
        p = ImmutablePolynomial(q,:y)
        BivariatePolynomial{F,X,N,Y,M}(p)
end
BivariatePolynomial(t::Tuple) = BivariatePolynomial{:x,:y}(t)
function BivariatePolynomial(r::ImmutablePolynomial{ImmutablePolynomial{F,X},Y,M}) where {F,X,Y,M}
    BivariatePolynomial(Tuple(z.coeffs for z in r.coeffs))
end

(p::BivariatePolynomial)(x,y) = p.p(y)(x)
(p::BivariatePolynomial)(x::AbstractVector) = p(x...)


for op in (:+,:-,:*)
    expr = Meta.parse("(Base.:$op)(p::BivariatePolynomial,q::BivariatePolynomial) = BivariatePolynomial($op(p.p,q.p))")
    eval(expr)
    expr = Meta.parse("(Base.:$op)(a::Number,q::BivariatePolynomial) = BivariatePolynomial($op(a,q.p))")
    eval(expr)
end


######################################
####      TensorPolynomial     ####
######################################
"""
```
    julia> ϕ = TensorPolynomial((2,1.,3.),(3.,2.))
    julia> ϕ(1,-1)
        6.0
    julia> ϕ([1,-1])
        6.0
```
Defines a scalar field `ϕ:ℝ²→ℝ` given by a product of a polynomial in `x` and  a polynomial in `y`. In this case: `(2.0 + 1.0x + 3.0x²)(3.0+2.0y)`.

The indeterminates can be specified:
```
    julia> ϕ = TensorPolynomial((2,1.,3.),(3.,2.),:z,:ξ)
```
produces the same polynomials, but on the variables `(z,ξ)`. 
"""
struct TensorPolynomial{F,X,N,Y,M} <: PolyScalarField{F,X,N,Y,M}
    px::ImmutablePolynomial{F,X,N}
    py::ImmutablePolynomial{F,Y,M}
end

function TensorPolynomial(t1::Tuple,t2::Tuple,X,Y)
    t1 = promote(t1...)
    t2 = promote(t2...)
    F = promote_type(eltype(t1),eltype(t2))
    N = length(t1)
    M = length(t2)
    p1 = ImmutablePolynomial(NTuple{N,F}(convert.(F,t1)),X)
    p2 = ImmutablePolynomial(NTuple{M,F}(convert.(F,t2)),Y)
    TensorPolynomial{F,X,N,Y,M}(p1,p2)
end
TensorPolynomial(t1::Tuple,t2::Tuple) =  TensorPolynomial(t1,t2,:x,:y)
(s::TensorPolynomial)(x,y) = s.px(x)*s.py(y)
(s::TensorPolynomial)(x::T) where T<:AbstractVector = s.px(x[1])*s.py(x[2])



Base.:*(a::Number,p::TensorPolynomial) = TensorPolynomial(a*p.px,p.py)

function Base.:*(p::AbstractPolynomial,q::TensorPolynomial)
    if indeterminate(p)===indeterminate(q.px)
        TensorPolynomial(p*q.px,q.py)
    elseif indeterminate(p)===indeterminate(q.py)
        TensorPolynomial(q.px,p*q.py)
    else
        throw(ArgumentError("Indeterminates does not match. Polynomial p has indeterminate $(indeterminate(p)), whereas TensorPolynomial q has indeterminates $(indeterminate(q.px)),$(indeterminate(q.py))."))
    end
end
Base.:*(q::TensorPolynomial,p::AbstractPolynomial) = p*q

function Base.:*(p::TensorPolynomial{F,X,N,Y,M},q::TensorPolynomial{F,X,K,Y,L}) where {F,X,N,Y,M,K,L}
    TensorPolynomial(p.px*q.px,p.py*q.py)
end


function Base.convert(BivariatePolynomial,p::TensorPolynomial{F,X,N,Y,M}) where {F,X,N,Y,M}
    cx = [p.px.coeffs...]
    cy = p.py.coeffs
    t = Tuple(Tuple(c*cx) for c in cy)
    BivariatePolynomial(t)
end
BivariatePolynomial(p::TensorPolynomial) = convert(BivariatePolynomial,p)
Base.promote(p::BivariatePolynomial,q::TensorPolynomial) = (p,BivariatePolynomial(q))
Base.promote(q::TensorPolynomial,p::BivariatePolynomial) = (BivariatePolynomial(q),p)


for op in (:+,:-,:*)
    expr = Meta.parse("(Base.:$op)(p::BivariatePolynomial,q::TensorPolynomial) = $op(p,convert(BivariatePolynomial,q))")
    eval(expr)
    expr = Meta.parse("(Base.:$op)(q::TensorPolynomial,p::BivariatePolynomial) = $op(p,q)")
    eval(expr)
end
for op in (:+,:-)
    expr = Meta.parse("(Base.:$op)(p::TensorPolynomial,q::TensorPolynomial) = $op(convert(BivariatePolynomial,p),convert(BivariatePolynomial,q))")
    eval(expr)
end


###############################
#       VECTOR FIELDS
###############################
struct PolyVectorField{F,X,N1,N2,Y,M1,M2} <: PolyField{F,X,Y}
    s1::PolyScalarField{F,X,N1,Y,M1}
    s2::PolyScalarField{F,X,N2,Y,M2}
end

(v::PolyVectorField)(x,y) = HPPoint(v.s1(x,y),v.s2(x,y))
(v::PolyVectorField)(x) = v(x[1],x[2])

Base.IteratorSize(::PolyVectorField) = HasLength()
Base.length(::PolyVectorField) = 2
get_x(p::PolyVectorField) = p.s1
get_y(p::PolyVectorField) = p.s2
function Base.getindex(p::PolyVectorField,i)
    if i == 1
        return get_x(p)
    elseif i == 2
        return get_y(p)
    else
        throw(ArgumentError("Index out of range. PolyVectorField has only 2 components."))
    end
end

function Base.:*(p::PolyScalarField,v::PolyVectorField)
    issubset(indeterminates(p),indeterminates(v)) || throw(ArgumentError("Fields have different indeterminates"))
    PolyVectorField(p*v[1],p*v[2])
end
function Base.:*(p::AbstractPolynomial,v::PolyVectorField)
    indeterminate(p) in indeterminates(v) || throw(ArgumentError("Indeterminates does not match."))
    PolyVectorField(p*v[1],p*v[2])    
end
Base.:*(a::Number,v::PolyVectorField) = PolyVectorField(a*v[1],a*v[2])
function Base.:*(w::AbstractMatrix,v::PolyVectorField)
    size(w) == (2,1) && return _vector_polyvector(w,v)
    size(w) == (2,2) && return _matrix_polyvector(w,v)
    throw(ArgumentError("It is only possible to multiply by a matrix of 2×2 or 2⨯1. For vector-vector field multiplication, check ⋅."))
end

_vector_polyvector(w,v) = w[1]*v[1] + w[2]*v[2]
_matrix_polyvector(w,v) = PolyVectorField(w[1,1]*v[1]+w[1,2]*v[2],w[2,1]*v[1]+w[2,2]*v[2])

function LinearAlgebra.dot(w::AbstractVector,v::PolyVectorField)
    length(w)!=2 && throw(ArgumentError("Vector must have length 2."))
    _vector_polyvector(w,v)
end
LinearAlgebra.dot(v::PolyVectorField,w::AbstractVector) = w⋅v

function Base.:*(w::AbstractVector,p::PolyField)
     length(w)!=2 && throw(ArgumentError("Vector must have length 2."))
     PolyVectorField(w[1]*p,w[2]*p)
end
Base.:*(p::PolyField,w::AbstractVector) = w*p

