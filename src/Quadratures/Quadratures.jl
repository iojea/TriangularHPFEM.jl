"""
    This module is a fork of `GrundmannMoeller.jl` (https://github.com/eschnett/GrundmannMoeller.jl.git). The only relevant difference is that we only consider 1D, 2D or 3D quadratures up to order 35, for which we have precomputed the size of the quadrature (number of points and weights), which allows us to store the data directly into StaticArrays, reducing the memory allocation. 
"""
module Quadratures 

using LinearAlgebra
using StaticArrays
using Dictionaries
using HyperDualNumbers
using TriangularhpFEM.Meshes

include("auxiliarydata.jl")
include("refelement.jl")
include("gmquads.jl")
include("measures.jl")

export QScheme, Measure, AuxDegData
export gmquadrature


end; #module

