module Spaces

using StaticArrays
using Dictionaries
using LinearAlgebra
using Polynomials
using ..Meshes
using ..PolyFields
#import ..PolyFields: gradient, divergence, laplacian
#include("basis.jl")
#include("poly.jl")
#include("basis.jl")
#include("opfields.jl")
include("hpspaces.jl")
#include("differentiation.jl")
include("operation.jl")



#export BasisIterator,StdScalarSpace,GradientSpace
# export TensorPolynomial
# export LegendreIterator
# export StandardBasis
export AbstractSpace
export StdScalarSpace, StdVectorSpace, StdTensorSpace
export ScalarSpace, VectorSpace, TensorSpace
export basis
export OperatorSpace
export gradient,divergence,laplacian
export dot
export Order
export order

end; #module
