module PolyFields

using Polynomials
using LinearAlgebra
using ..Meshes

include("fields.jl")
include("opfields.jl")
include("legendre.jl")
include("differentiation.jl")
include("integration.jl")
include("show.jl")

const Δ = laplacian    
const ∇ = gradient

export BivariatePolynomial,TensorPolynomial
export PolyScalarField, PolyVectorField
export FieldType
export OperationField
export LegendreIterator,StandardBasis
export derivative,divergence,laplacian,gradient,∇,Δ
export ref_integrate
export printpoly2,show


end; #module
