module Meshes

using Dictionaries
using ExactPredicates
using LinearAlgebra
using Makie
using Markdown
using Printf
using Triangulate
using StaticArrays
using DocStringExtensions


export HPTriangulation,HPPoint,HPEdge,HPTriangle,HPMesh
export meshhp,circmesh,circmesh_graded_center,rectmesh
export show, plothpmesh, animate_refinement


include("tuplehp.jl")
include("point.jl")
include("edge.jl")
include("triangle.jl")
include("dofs.jl")
include("mesh.jl")
include("refine.jl")
include("show.jl")
include("plots.jl")
include("examples.jl")

const HASH_SEED = UInt === UInt64 ? 0x793bac59abf9a1da : 0xdea7f1da

end;