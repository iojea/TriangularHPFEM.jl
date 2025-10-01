
abstract type FieldsType end
struct PolyType <: FieldsType end
struct GeneralType <:FieldsType end

struct OperationField{F<:Function,F1,F2}
    op::F
    args::Tuple{F1,F2}
end

arguments(of::OperationField) = of.args
operation(of::OperationField) = of.op

(of::OperationField)(x) = operation(of)((f(x) for f in arguments(of))...)

FieldsType(::PolyField) = PolyType()
Base.promote_type(::PolyType,::PolyType) = PolyType()
Base.promote_type(::PolyType,::GeneralType) = GeneralType()
Base.promote_type(::GeneralType,::PolyType) = GeneralType()
Base.promote_type(::GeneralType,::GeneralType) = GeneralType()
FieldsType(::A,::B) where {A,B} = promote_type(A,B)
FieldsType(op::OperationField) = promote_type(FieldsType.(op.args)...)

Base.:+(p::P,q::P) where {P<:PolyField} = OperationField(+,(p,q))
Base.:-(p::P,q::P) where {P<:PolyField} = OperationField(-,(p,q))

function Base.:*(A::T,u::PolyVectorField) where {T<:AbstractMatrix}
    println(size(A))
    println(typeof(A))
    size(A) == (2,2) && return matpolyprod(A,u)
    size(A) == (2,1) && return vecpolyprod(A[:],u)
    throw(ArgumentError("Sizes does not match"))
end

function matpolyprod(A::T,u::PolyVectorField) where {T<:AbstractMatrix}
    o1 = OperationField(+,(A[1,1]*u.s1,A[1,2]*u.s2))
    o2 = OperationField(+,(A[2,1]*u.s1,A[2,2]*u.s2))
    OperationField(+,(o1,o2))
end

function vecpolyprod(v::AbstractVector,u::PolyVectorField)
    length(v)== 2 || throw(ArgumentError("Sizes does not match."))
    OperationField(+,(v[1]*u.s1,v[2]*u.s2))
end

LinearAlgebra.:⋅(v::AbstractVector,u::PolyVectorField) = vecpolyprod(v,u)
LinearAlgebra.:⋅(A::AbstractMatrix,u::PolyVectorField) = *(A,u)
function LinearAlgebra.:⋅(u1::PolyVectorField,u2::PolyVectorField)
    OperationField(+,(u1.s1*u2.s1,u1.s2*u2.s2))
end


Base.:*(f::Function,u::PolyScalarField) = OperationField(*,(f,u))
Base.:*(u::PolyScalarField,f::Function) = f*u





