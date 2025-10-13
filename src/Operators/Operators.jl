module Operators

using LinearAlgebra
using Polynomials
using ..Spaces
using ..PolyFields

include("differentiation.jl")
include("operation.jl")
const Δ = laplacian    
const ∇ = gradient

export gradient, divergence,∇
export Operation
export dot
end; #module