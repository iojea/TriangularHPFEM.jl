module PolyFields

using Polynomials
using LinearAlgebra
using ..Meshes

include("fields.jl")
include("opfields.jl")
include("legendre.jl")
include("differentiation.jl")
include("show.jl")

const Δ = laplacian    
const ∇ = gradient

export BivariatePolynomial,TensorPolynomial
export PolyScalarField, PolyVectorField
export OperationField
export LegendreIterator,StandardBasis
export derivative,divergence,laplacian,gradient,∇,Δ


end; #module
