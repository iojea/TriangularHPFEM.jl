abstract type PolyField{F,X,Y} <:Function end


###########################################################################
#################              SCALAR FIELDS              #################
###########################################################################
#
abstract type PolyScalarField{F,X,Y} <: PolyField{F,X,Y} end
indeterminate(::AbstractPolynomial{T,X}) where {T,X} = X
indeterminates(::PolyField{F,X,Y}) where {F,X,Y} = (X,Y)
Base.length(::PolyScalarField) = 1
# 

_zerofill(t::NTuple{K,T},N) where {K,T} = K>=N ? t : (t...,zeros(T,N-K)...)
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
struct TensorPolynomial{F,X,Y} <: PolyScalarField{F,X,Y}
    px::ImmutablePolynomial{F,X,N} where N
    py::ImmutablePolynomial{F,Y,M} where M
end

function TensorPolynomial(t1::Tuple,t2::Tuple,X,Y)
    t1 = promote(t1...)
    t2 = promote(t2...)
    F = promote_type(eltype(t1),eltype(t2))
    N = length(t1)
    M = length(t2)
    p1 = ImmutablePolynomial(NTuple{N,F}(convert.(F,t1)),X)
    p2 = ImmutablePolynomial(NTuple{M,F}(convert.(F,t2)),Y)
    TensorPolynomial{F,X,Y}(p1,p2)
end
TensorPolynomial(t1::Tuple,t2::Tuple) =  TensorPolynomial(t1,t2,:x,:y)

(s::TensorPolynomial)(x,y) = s.px(x)*s.py(y)
(s::TensorPolynomial)(x::T) where T<:AbstractVector = s.px(x[1])*s.py(x[2])


Base.promote(t::TensorPolynomial{F,X,Y},s::TensorPolynomial{F,X,Y}) where {F,X,Y}= (t,s)
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

function Base.:*(p::TensorPolynomial{F,X,Y},q::TensorPolynomial{F,X,Y}) where {F,X,Y}
    TensorPolynomial(p.px*q.px,p.py*q.py)
end
LinearAlgebra.:⋅(p::TensorPolynomial,q::TensorPolynomial) = p*q
function Base.zero(::TensorPolynomial{F,X,Y}) where {F,X,Y}
    TensorPolynomial((zero(F),),(zero(F),))
end
function Base.zero(::Type{TensorPolynomial{F,X,Y}}) where {F,X,Y}
    TensorPolynomial((zero(F),),(zero(F),))
end

###############################
#       TENSOR FIELDS
###############################
struct PolyTensorField{F,X,Y,T<:PolyScalarField{F,X,Y},N} <: PolyField{F,X,Y}
    tensor::Array{T,N}
end
function PolyTensorField{T,N}() where {N,T<:PolyScalarField}
    a = Array{T,N}(undef,[2 for _ in 1:N]...)
    PolyTensorField(a)
end

# PolyVectorField
const PolyVectorField{F,X,Y,T} = PolyTensorField{F,X,Y,T,1}
PolyVectorField(x::AbstractArray{T,1}) where T = PolyTensorField(x)

#PolyMatrixField
const PolyMatrixField{F,X,Y,T} = PolyTensorField{F,X,Y,T,2}
PolyMatrixField(x::AbstractArray{T,2}) where T = PolyTensorField(x)
function (v::PolyTensorField{F,X,Y,T,N})(x,y) where {F,X,Y,T,N}
    z = Array{F,N}(undef,[2 for _ in 1:N]...)
    for i in eachindex(v.tensor)
        z[i] = v.tensor[i](x,y)
    end
    z
end


function (v::PolyTensorField{F,X,Y,T,N})(x) where {F,X,Y,T,N}
    v(x[1],x[2])
end

Base.IteratorSize(::PolyTensorField{F,X,Y,T,N}) where {F,X,Y,T,N} = Base.HasShape{N}()
Base.length(p::PolyTensorField) = length(p.tensor)
Base.size(p::PolyTensorField) = size(p.tensor)
Base.iterate(p::PolyTensorField,st=nothing) = iterate(p.tensor,st)
Base.getindex(p::PolyTensorField,i) = getindex(p.tensor,i)

Base.:*(a::Number,p::PolyTensorField) = PolyTensorField(a*p.tensor)
Base.:*(p::PolyTensorField,a::Number) = a*p

# This function needs improvement. c is created with an AbstractType but it should be create with ConcreteType when possible. 
function Base.:*(A::AbstractMatrix,p::PolyTensorField{F,X,Y,T,N}) where {F,X,Y,T,N}
    sA = size(A); sp = size(p)
    sA[end] == sp[1] || throw(DimensionMismatch())
    D = length(sA)+length(sp)-2
    c = Array{PolyScalarField{F,X,Y},D}(undef,[2 for _ in 1:D]...)
    c .= A*p.tensor
    length(c) == 1 ? c : PolyTensorField(c)
end

function Base.:*(p::PolyScalarField,v::T) where T<:PolyTensorField
    issubset(indeterminates(p),indeterminates(v)) || throw(ArgumentError("Fields have different indeterminates"))
    T(p .* v.tensor)
end
Base.:*(v::T,p::PolyScalarField) where T<:PolyTensorField = p*v


function _outer(p::PolyVectorField,q::PolyVectorField)
    @cast a[i,j] := p.tensor[i]*q.tensor[j]
    PolyTensorField(a)
end

LinearAlgebra.dot(p::PolyVectorField,q::PolyVectorField) = dot(p.tensor,q.tensor)


struct Field
    op
    args::Tuple
end

Base.:*(f::Function,p::PolyField) = Field(*,(f,p))

