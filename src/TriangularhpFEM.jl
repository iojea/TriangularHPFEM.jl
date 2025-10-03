module TriangularhpFEM

#using DocStrigExtensions

include("Meshes/Meshes.jl")
#include("Quadratures/Quadratures.jl")
#include("Spaces/Spaces.jl")
include("PolyFields/PolyFields.jl")

macro publish(mod,name)
  quote
    using TriangularhpFEM.$mod: $name; export $name
  end
end

@publish Meshes HPPoint
@publish Meshes HPEdge
@publish Meshes HPTriangle
@publish Meshes HPMesh
@publish Meshes meshhp
@publish Meshes circmesh
@publish Meshes circmesh_graded_center
@publish Meshes rectmesh
@publish Meshes show
@publish Meshes plothpmesh
@publish Meshes animate_refinement

@publish PolyFields BivariatePolynomial
@publish PolyFields TensorPolynomial
@publish PolyFields PolyScalarField
@publish PolyFields PolyVectorField
@publish PolyFields OperationField
@publish PolyFields gradient
@publish PolyFields divergence
@publish PolyFields ∇
@publish PolyFields dot
@publish PolyFields laplacian
@publish PolyFields ref_integrate
@publish PolyFields LegendreIterator
@publish PolyFields StandardBasis
end;

#@publish Quadratures QScheme
#@publish Quadratures gmquadrature
#@publish Quadratures hatmass
#@publish Quadratures hatstiff
#@publish Quadratures compute_dimension
#@publish Quadratures boundary_nodes
#@publish Quadratures matrix_F
#@publish Quadratures matrix_C
#@publish Quadratures degs


#@publish Spaces BasisIterator
#@publish Spaces StdScalarSpace
#@publish Spaces GradientSpace
#@publish Spaces Iterator
#@publish Spaces Integrand
#@publish Spaces ∫
#@publish Spaces gradient
#@publish Spaces ∇
#@publish Spaces get_space
#@publish Spaces get_spaces

#end; #module
