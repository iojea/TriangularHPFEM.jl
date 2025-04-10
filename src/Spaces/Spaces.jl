module Spaces

using StaticArrays
using Dictionaries
using LinearAlgebra
using TriangularhpFEM.Meshes
using TriangularhpFEM.Quadratures


include("basis.jl")
include("dofs.jl")
include("hpspaces.jl")

export BasisIterator,StdScalarSpace,GradientSpace
export Integrand
export gradient,∇
export get_space,get_spaces

end; #module
