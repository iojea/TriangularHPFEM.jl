module TriangularhpFEM

#using DocStrigExtensions

macro publish(mod,name)
  quote
    using TriangularhpFEM.$mod: $name; export $name
  end
end


include("Meshes/Meshes.jl")
# include("PolyFields/PolyFields.jl")


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
@publish Meshes circmesh
@publish Meshes circmesh_graded_center
@publish Meshes rectmesh
@publish Meshes show
@publish Meshes plothpmesh
@publish Meshes animate_refinement


@publish Quadratures QScheme
@publish Quadratures Measure
@publish Quadratures AuxDegData
@publish Quadratures AuxMeshData
@publish Quadratures gmquadrature
@publish Quadratures compute_dimension
@publish Quadratures boundartrilisty_nodes
@publish Quadratures matrix_F
@publish Quadratures matrix_C
@publish Quadratures degs
@publish Quadratures integrate
@publish Quadratures ref_integrate

@publish Spaces BivariatePolynomial
@publish Spaces TensorPolynomial
@publish Spaces PolyVectorField
@publish Spaces PolySum
@publish Spaces LegendreIterator
@publish Spaces StandardBasis
@publish Spaces PolySum
@publish Spaces show
@publish Spaces gradient
@publish Spaces divergence
@publish Spaces laplacian
@publish Spaces ∇

@publish Spaces Integrand
@publish Spaces StdScalarSpace
@publish Spaces StdVectorSpace
@publish Spaces StdTensorSpace
@publish Spaces ScalarSpace
@publish Spaces VectorSpace
@publish Spaces TensorSpace
@publish Spaces OperatorSpace
@publish Spaces Operation
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
