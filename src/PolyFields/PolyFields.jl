module PolyFields

using Polynomials
using LinearAlgebra
using TensorCast
using ..Meshes

include("fields.jl")
include("opfields.jl")
include("differentiation.jl")
include("legendre.jl")

include("show.jl")

const ∇ = gradient
const Δ = laplacian
const ⊗ = _outer
export PolyField
export TensorPolynomial
export PolyScalarField, PolyVectorField, PolyTensorField
export Field
export PolySum
export LegendreIterator,StandardBasis
export gradient, divergence, laplacian,_outer
export dot

end; #module
