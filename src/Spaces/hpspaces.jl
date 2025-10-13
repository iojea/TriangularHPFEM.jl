
abstract type AbstractSpace end

abstract type ScalarSpace <: AbstractSpace end
abstract type VectorSpace <: AbstractSpace end
abstract type TensorSpace <: AbstractSpace end

struct StdScalarSpace <: ScalarSpace end
struct StdVectorSpace <: VectorSpace end
struct StdTensorSpace <: TensorSpace end

struct OperatorSpace{F<:Function,S<:AbstractSpace} <: AbstractSpace
    operator::F
    space::S
end



# Traits for
abstract type SpaceType end
struct ScalarSpace <: SpaceType end
struct VectorSpace <: SpaceType end
struct TensorSpace <: SpaceType end

spacetype(s::StdScalarSpace) = ScalarSpace
spacetype(s::StdVectorSpace) = VectorSpace
spacetype(s::StdTensorSpace) = TensorSpace

gradient(s::StdScalarSpace) = OperatorSpace(gradient,s)
divergence(s::StdVectorSpace) = OperatorSpace(divergence,s)
laplacian(s::StdScalarSpace) = OperatorSpace(laplacian,s)


#spacetype(s::OperatorSpace{::Type{typeof(gradient),StdScalarSpace}}) = VectorSpace

# LinearAlgebra.dot()

# Base.:*(a::ScalarHPSpace,b::AbstractHPSpace) = OperationField(*,(a,b))
# Base.:*(m::Array,a::VectorHPSpace) = OperationField(*,(m,a))
# Base.:*(f::Function,a::VectorHPSpace) = OperationField(*,(f,a))





# A trait for identifying Constant Coefficients, which allow precomputation of local tensors.
# abstract type CoeffType end

# struct Constant <: CoeffType end
# struct Variable <: CoeffType end

# coefftype(::AbstractHPSpace) = Constant()
# coefftype(::AbstractArray) = Constant()
# coefftype(::Number) = Constant()

# coefftype(::Function) = Variable()
# coefftype(p::PolyField) = degree(p)==0 ? Constant() : Variable()

# # Sintactic sugar
# struct Integrand{F<:Function,F1,F2}
#     op::F
#     args::Tuple{F1,F2}
# end
# const ∫ = Integrand

# Base.:*(f::Integrand,m::Measure) = integrate(CoeffType(f.func),f.func,m)


# get_space(space::AbstractHPSpace) = space
# function get_space(opfield::OperationField)
#     (;args) = opfield
#     if typeof(args[1])<:AbstractHPSpace
#         return args[1]
#     elseif typeof(args[2])<:AbstractHPSpace
#         return args[2]
#     else
#         error("Something is wrong. An OperatorField was found with no space involved")
#     end
# end
# get_spaces(itd::Integrand) = get_space.(itd.func.args)

# function integrate(::Constant,f::Integrand,m::Measure)
#     M = build_tensors(f,m)
# end


# #function build_tensors()