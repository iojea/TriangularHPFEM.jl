module TriangularhpFEM

#using DocStrigExtensions

macro publish(mod,name)
  quote
    using TriangularhpFEM.$mod: $name; export $name
  end
end


include("Meshes/Meshes.jl")
include("PolyFields/PolyFields.jl")
include("Spaces/Spaces.jl")
include("Quadratures/Quadratures.jl")



@publish Meshes HPPoint
@publish Meshes HPEdge
@publish Meshes HPTriangle
@publish Meshes HPMesh
@publish Meshes meshhp
@publish Meshes triangles
@publish Meshes edges
@publish Meshes degrees
@publish Meshes degree
@publish Meshes degrees_of_freedom!
@publish Meshes circmesh
@publish Meshes circmesh_graded_center
@publish Meshes rectmesh
@publish Meshes show
@publish Meshes plothpmesh
@publish Meshes animate_refinement

@publish PolyFields PolyField
@publish PolyFields TensorPolynomial
@publish PolyFields PolyScalarField
@publish PolyFields PolyVectorField
@publish PolyFields PolySum
@publish PolyFields PolyTensorField 
@publish PolyFields gradient
@publish PolyFields divergence
@publish PolyFields laplacian
@publish PolyFields ∇
@publish PolyFields LegendreIterator
@publish PolyFields StandardBasis
@publish PolyFields _outer
@publish PolyFields dot

@publish Quadratures QScheme
@publish Quadratures Measure
@publish Quadratures AuxDegData
@publish Quadratures AuxMeshData
@publish Quadratures gmquadrature
@publish Quadratures compute_dimension
@publish Quadratures boundary_nodes
@publish Quadratures matrix_F
@publish Quadratures matrix_C
@publish Quadratures degs
@publish Quadratures integrate
@publish Quadratures ref_integrate

@publish Spaces show
@publish Spaces polynize
@publish Spaces basis

@publish Spaces Integrand
@publish Spaces StdScalarSpace
@publish Spaces StdVectorSpace
@publish Spaces StdTensorSpace
@publish Spaces ScalarSpace
@publish Spaces VectorSpace
@publish Spaces TensorSpace
@publish Spaces OperatorSpace
@publish Spaces Operation
@publish Spaces gradient
@publish Spaces divergence
@publish Spaces laplacian
@publish Spaces ∇
@publish Spaces coefftype
@publish Spaces dot
@publish Spaces ∫
@publish Spaces order
end;

#@publish Spaces BasisIterator
#@publish Spaces StdScalarSpace
#@publish Spaces GradientSpace
#@publish Spaces Iterator
#@publish Spaces Integrand
#@publish Spaces gradient
#@publish Spaces ∇
#@publish Spaces get_space
#@publish Spaces get_spaces

#end; #module
