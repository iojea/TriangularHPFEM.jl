module PolyFields

using Polynomials
using LinearAlgebra
using ..Meshes

include("fields.jl")
include("opfields.jl")
include("differentiation.jl")
include("legendre.jl")

include("show.jl")


export BivariatePolynomial,TensorPolynomial
export PolyScalarField, PolyVectorField
export OperationField
export LegendreIterator,StandardBasis
export gradient, divergence, laplacian

end; #module
