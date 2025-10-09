abstract type AbstractHPSpace end

abstract type ScalarHPSpace <: AbstractHPSpace end

abstract type VectorHPSpace <: AbstractHPSpace end

abstract type TensorHPSpace <: AbstractHPSpace end

struct StdScalarSpace <: ScalarHPSpace end


struct TestSpace{T<:AbstractHPSpace} <: AbstractHPSpace
    space::T
end

struct TrialSpace{T<:AbstractHPSpace,F<:Function} <:AbstractHPSpace
    space::T
    g::F
end


struct OperatorSpace{F<:Function,S<:AbstractHPSpace}
    op::F
    space::S
end

LinearAlgebra.:⋅(a::VectorHPSpace,b::VectorHPSpace) = OperationField(⋅,(a,b))
Base.:*(a::ScalarHPSpace,b::AbstractHPSpace) = OperationField(*,(a,b))
Base.:*(m::Array,a::VectorHPSpace) = OperationField(*,(m,a))Base.:*(f::Function,a::VectorHPSpace) = OperationField(*,(f,a))





# A trait for identifying Constant Coefficients, which allow precomputation of local tensors.
abstract type CoeffType end

struct Constant <: CoeffType end
struct Variable <: CoeffType end

coefftype(::AbstractHPSpace) = Constant()
coefftype(::AbstractArray) = Constant()
coefftype(::Number) = Constant()

coefftype(::Function) = Variable()
coefftype(p::PolyField) = degree(p)==0 ? Constant() : Variable()

# Sintactic sugar
struct Integrand{F<:Function,F1,F2}
    op::F
    args::Tuple{F1,F2}
end
const ∫ = Integrand

Base.:*(f::Integrand,m::Measure) = integrate(CoeffType(f.func),f.func,m)


get_space(space::AbstractHPSpace) = space
function get_space(opfield::OperationField)
    (;args) = opfield
    if typeof(args[1])<:AbstractHPSpace
        return args[1]
    elseif typeof(args[2])<:AbstractHPSpace
        return args[2]
    else
        error("Something is wrong. An OperatorField was found with no space involved")
    end
end
get_spaces(itd::Integrand) = get_space.(itd.func.args)

function integrate(::Constant,f::Integrand,m::Measure)
    M = build_tensors(f,m)
end


#function build_tensors()