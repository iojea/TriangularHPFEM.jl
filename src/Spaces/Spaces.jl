module Spaces

using StaticArrays
using Dictionaries
using LinearAlgebra
using ..Meshes
using ..Quadratures
using ..PolyFields
import ..PolyFields: gradient, divergence, laplacian
#include("basis.jl")
include("hpspaces.jl")
include("operation.jl")

#export BasisIterator,StdScalarSpace,GradientSpace
export Integrand
export AbstractSpace
export StdScalarSpace, StdVectorSpace, StdTensorSpace
export ScalarSpace, VectorSpace, TensorSpace
export OperatorSpace
export spacetype
export gradient,divergence,laplacian

end; #module
