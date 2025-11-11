module Matrices

using StaticArrays
using LinearAlgebra
using ..Meshes
using ..PolyFields
using ..Spaces
using ..Quadratures


struct AffineTransformation{F}
    A::MMatrix{F,2,2}
    iA::MMatrix{F,2,2}
    b::MVector{F,2}
    jacobian::Base.RefValue{F}
    function AffineTransform{F}(A,b) where F<:Number
        @assert size(A)==(2,2)
        @assert size(b)==(2,)
        new{F}(F.(A),F.(inv(A)),F.(b),Base.RefValue(abs(det(A))))
    end
end
function AffineTransform(A,b)
    TA = eltype(A)
    Tb = eltype(b)
    T = promote_type(TA,Tb)
    AffineTransform{T}(A,b)
end
AffineTransform{F}() where F = AffineTransform(@MMatrix(zeros(F,2,2)),@MVector(zeros(F,2)))

(aff::AffineTransform)(x) = aff.A*x+aff.b

function update!(aff::AffineTransform,vert)
    (;A,iA,b,jacobian) = aff
    @views A[:, 1] .= 0.5(vert[:, 3] - vert[:, 2])
    @views A[:, 2] .= 0.5(vert[:, 1] - vert[:, 2])
    @views b .= 0.5(vert[:, 1] + vert[:, 3])
    iA .= inv(A)
    jacobian[] = abs(det(A))
end

# A trait for evaluation of Field
abstract type EvalType end
struct Eval <: EvalType end
struct Compose <: EvalType end
struct Pass <: EvalType end

evaltype(_) = Pass()
evaltype(::Function) = Compose()
evaltype(::PolyField) = Eval()

evaluate(::Eval,f,x,t::AffineTransform) = f(x)
evaluate(::Field,f,x,t::AffineTransform) = f(t(x))
evaluate(::Pass,f,x,t) = f

(o::Field)(x,t) = o.op((evaluate(evaltype(arg),arg,x,t) for arg in o.args)...)
end;