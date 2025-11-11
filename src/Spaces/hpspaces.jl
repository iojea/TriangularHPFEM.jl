
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


PolyFields.gradient(s::StdScalarSpace) = OperatorSpace(gradient,s)
PolyFields.divergence(s::StdVectorSpace) = OperatorSpace(divergence,s)
PolyFields.laplacian(s::StdScalarSpace) = OperatorSpace(laplacian,s)


# A trait for identifying Constant Coefficients, which allow precomputation of local tensors.
abstract type CoeffType end

struct Constant <: CoeffType end
struct Variable <: CoeffType end

coefftype(::AbstractSpace) = Constant
coefftype(::AbstractArray) = Constant
coefftype(::Number) = Constant

coefftype(::Function) = Variable

Base.promote_type(::Union{Constant,Type{Constant}},::Union{Variable,Type{Variable}}) = Variable
Base.promote_type(::Union{Constant,Type{Constant}},::Union{Constant,Type{Constant}}) = Constant
Base.promote_type(::Union{Variable,Type{Variable}},::Union{Constant,Type{Constant}}) = Variable
Base.promote_type(::Union{Constant,Type{Constant}},::Union{Nothing,Type{Nothing}}) = Constant
Base.promote_type(::Union{Nothing,Type{Nothing}},::Union{Constant,Type{Constant}}) = Constant
Base.promote_type(::Union{Nothing,Type{Nothing}},::Union{Variable,Type{Variable}}) = Variable
Base.promote_type(::Union{Variable,Type{Variable}},::Union{Nothing,Type{Nothing}}) = Variable


struct Order{B} end

order(::Number) = Order{0}()
order(::AbstractArray) = Order{0}()
order(::PolyField) = Order{0}()
order(::Function) = Order{0}()
order(::AbstractSpace) = Order{(0,)}()
order(::OperatorSpace{typeof(gradient),S}) where S = Order{(1,)}()
order(::OperatorSpace{typeof(divergence),S}) where S = Order{(1,)}()
order(::OperatorSpace{typeof(laplacian),S}) where S = Order{(2,)}()
combine(::Order{B},::Order{C}) where {B,C} = Order{(B...,C...)}()
combine(::Order{B},::Order{0}) where B = Order{B}()
combine(::Order{0},::Order{B}) where B = Order{B}()


basis(::StdScalarSpace) = StandardBasis

polynize(::AbstractSpace,p::PolyField) = p
polynize(op::OperatorSpace,p::PolyField) = op.operator(p)

# # Sintactic sugar

struct Integrand
    op
end
#Integrand(op) = Integrand()

# Base.:*(f::Integrand,m::Measure) = integrate(CoeffType(f.op),f.op,m)


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