abstract type AbstractHPSpace end

abstract type ScalarHPSpace <: AbstractHPSpace end

abstract type VectorHPSpace <: AbstractHPSpace end

abstract type TensorHPSpace <: AbstractHPSpace end

struct StdScalarSpace <: ScalarHPSpace end

struct GradientSpace <: VectorHPSpace end

gradient(::StdScalarSpace) = GradientSpace()

const ∇ = gradient


struct TestSpace{T<:AbstractHPSpace} <: AbstractHPSpace
    space::T
    dofs::DOFs
end

struct TrialSpace{T<:AbstractHPSpace,F<:Function} <:AbstractHPSpace
    space::T
    g::F
end


struct OperationField{F<:Function,T<:Tuple}
    op::F
    args::T
end

LinearAlgebra.:⋅(a::VectorHPSpace,b::VectorHPSpace) = OperationField(⋅,(a,b))
Base.:*(a::ScalarHPSpace,b::AbstractHPSpace) = OperationField(*,(a,b))



# A trait for identifying Constant Coefficients, which allow precomputation of local tensors.
abstract type CoeffType end

struct Constant <: CoeffType end
struct Variable <: CoeffType end

# Sintactic sugar
struct Integrand
    func
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