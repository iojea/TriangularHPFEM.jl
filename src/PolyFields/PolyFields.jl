module PolyFields

using Polynomials
using LinearAlgebra
using TensorOperations
using ..Meshes

include("fields.jl")
include("opfields.jl")
include("differentiation.jl")
include("legendre.jl")

include("show.jl")

const ∇ = gradient
const Δ = laplacian
export PolyField
export TensorPolynomial
export PolyScalarField, PolyVectorField, PolyTensorField
export PolySum
export LegendreIterator,StandardBasis
export gradient, divergence, laplacian,outer

end; #module
