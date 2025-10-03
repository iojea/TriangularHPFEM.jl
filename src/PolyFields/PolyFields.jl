module PolyFields

using Polynomials
using LinearAlgebra

include("poly2d.jl")
include("fields.jl")
include("opfields.jl")
include("legendre.jl")
include("differentiation.jl")
include("integration.jl")

const Δ = laplacian    
const ∇ = gradient

export BivariatePolynomial,TensorPolynomial
export PolyScalarField, PolyVectorField
export OperationField
export LegendreIterator,StandardBasis
export derivative,divergence,laplacian,gradient,∇,Δ
export ref_integrate


end; #module
