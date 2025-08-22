const TriangleList{I,P,F} = Dictionary{HPTriangle{I},TriangleProperties{P,F}} where {I,P,F}
const EdgeList{I,P} = Dictionary{HPEdge{I},EdgeProperties{P,Bool}} where {I,P}
#const AuxList{P,N,F} = Dictionary{NTuple{3,P},AuxData{N,F}} where {P,N,F}


abstract type HPTriangulation end


"""
    $(SIGNATURES)

A mesh for `HP` finite element methods. Its fields are: 
    + `points::Vector{HPPoint{2,F}}`: a vector of points
    + `trilist::TriangleList{I,P,F}`: set of trilist
    + `edgelist::EdgeList{I,P}`: set of edgelist
    + `auxlist::AuxList{P,N,F}`: auxiliary data for integrating the 

    HPMesh(tri::TriangulationIO)
builds an `HPMesh` from a `Triangulate.TriangulatioIO` struct.

Meshes of type `HPMesh` can also be constructed using the helper functions such us `circmesh`, `rectmesh`, `squaremesh`. For a more general constructor, check the docs for `meshhp`.
"""
struct HPMesh{F<:AbstractFloat,I<:Integer,P<:Integer} <: HPTriangulation
    points::Vector{HPPoint{2,F}} #ElasticMatrix{F,Vector{F}}
    trilist::TriangleList{I,P,F}
    edgelist::EdgeList{I,P}
    #auxdata::AuxList{P,N,F}   
end

function HPMesh(mat::Matrix,tris::TriangleList,edgs::EdgeList)
    HPMesh(HPPoint.(eachcol(mat)),tris,edgs)
end

function maybeconvert(::Type{T},arr::AbstractArray) where T
    eltype(arr)==T ? arr : convert.(T,arr)
end

function HPMesh{F,I,P}(tri::TriangulateIO) where {F,I,P}
    (;pointlist,trianglelist,edgelist,edgemarkerlist) = tri
    points = HPPoint{2,F}.(point for point in eachcol(pointlist))
    edgelist = maybeconvert(I,edgelist)
    triangles = dictionary([triangle(I,t,pointlist) => TriangleProperties(P,F) for t in eachcol(trianglelist)])
    edges = dictionary([HPEdge(e) => EdgeProperties(one(P),P(edgemarkerlist[i]),false) for (i,e) in enumerate(eachcol(edgelist))] )
    HPMesh(points,triangles,edges)
end

function Base.copy(mesh::HPMesh)
    HPMesh(deepcopy(mesh.points),deepcopy(mesh.trilist),deepcopy(mesh.edgelist))
end

@inline edges(list::T) where T<:EdgeList = keys(list)
@inline triangles(list::T) where T<:TriangleList = keys(list)
@inline edges(mesh::T) where T<:HPMesh = keys(mesh.edgelist)
@inline triangles(mesh::T) where T<:HPMesh = keys(mesh.trilist) 



@inline isgreen(t::HPTriangle,m::HPMesh) = m.trilist[t].refine[] == 1
@inline isblue(t::HPTriangle,m::HPMesh)  = m.trilist[t].refine[] == 2
@inline isred(t::HPTriangle,m::HPMesh)  = m.trilist[t].refine[] == 3

"""
    $(SIGNATURES)

Checks if 'e' is stored as presented or in reverse order. 
"""
function same_order(e::HPEdge{I},elist::EdgeList{I,P}) where {I,P}
    _,(_,k) = gettoken(elist,e)
    oe      = gettokenvalue(keys(elist),k)
    oe == e
end

"""
    $(SIGNATURES)

returns a permutation of contiguous indices `ind` such that `v[ind]` is ordered. 
"""
function psortperm(v)
    i    = argmin(v)
    if v[i] ≤ v[mod1(i+1,3)] ≤ v[mod1(i+2,3)]
        ind = mod1.(i:i+2,3)
    else
        ind = mod1.(i:-1:i-2,3)
    end
    return ind
end

"""
  $(SIGNATURES)

Given a triangle `t` belonging to a `mesh`, it returns the (sorted) degrees
`p₁<=p₂<=p₃` of its edges, and the edges also sorted according to the degrees.
"""
function pedges(t::HPTriangle{I},mesh::HPMesh{F,I,P}) where {F,I,P}
    p,nod = pnodes(t,mesh)
    eds   = (HPEdge(nod[SVector(i,mod1(i+1,3))]) for i in 1:3)
    p,eds
end

"""
   $(SIGNATURES)

Given a triangle `t` belonging to a `mesh`, it returns the degrees `p₁<=p₂<=p₃` of the edges (sorted) and the nodes of `t` also sorted according to the degrees. 
"""
function pnodes(t::HPTriangle{I},mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    eds = edges(t)
    p   = degree.(getindices(edgelist,eds))
    pind = psortperm(p)
    p[pind],first.(eds[pind])
end

"""
    $(SIGNATURES)

returns a list containing the indices of the vertices of `mesh` tha lie in its boundary. 
"""
function boundarynodes(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    v = zeros(I,sum(>(0),marker.(edgelist)))
    i = 1
    for e in edges(edgelist)
        if marker(edgelist[e])>0
            for j in e
                if j ∉ v
                    v[i] = j
                    i   += 1
                end
            end
        end
    end
    v
end

"""
    $(SIGNATURES)

A helper function for creating `hp`-meshes from boundary data.
For simple mesh, it is enough matrix of `2×N` with the vertices of the polygon and a mesh size `h`:

```julia
julia> vertices = [-1. 0.;1. 0.;1.5 1.;0. 1.5;-1. 1.]'
julia> m = meshhp(vertices,0.1)
```

For more complex meshes, a matrix of `segments` and a vector of `markers` are needed. `segments` is a matrix of integers with size `2×S` indicating how vertices should be joined. `markers` is a vector of integers that impose a mark on each segment. The primary goal of markers is to indicate if a piece of boundary will hold Dirichlet (odd marker) or Neumann (even marker) conditions. If ommited, Dirichlet conditions will be assumed. In the following example we create a mesh of a square where Neumann conditions are imposed on the upper half.

```julia
julia> vert = [0. 0.;1. 0.;1. 0.5;1. 1.;0. 1.;0. 0.5]'
julia> segs = [1 2;2 3;3 4;4 5;5 6;6 1]'
julia> mark = [1,1,2,2,2,1]
julia> m = meshhp(vert,0.1;segments=segs,markers=mark)
```

Furthermore, meshes with holes can also be constructed. In this case, the `segments` should include the segments corresponding to the interior boudaries (always in a positively oriented order), and a parameter `holes` should be passed. `holes` is a matrix where each column contains a point included in one hole. Only one point per hole is necessary. 

```julia
julia> vert = [0. 0.;1. 0.;1. 0.5;1. 1.;0. 1.;0. 0.5;0.25 0.25;0.5 0.25;0.5 0.5;0.25 0.5]'
julia> segs = [1 2;2 3;3 4;4 5;5 6;6 1;7 8;8 9;9 10;10 7]'
julia> mark = [1,1,2,2,2,1,1,1,2,2]
julia> hole = [0.3 0.3]'
julia> m = meshhp(vert,0.1;segments=segs,markers=mark,holes=hole)
```
"""
function meshhp(vertices,h;segments=nothing,markers=nothing,holes=nothing)
    tri = TriangulateIO()
    tri.pointlist = vertices
    if isnothing(segments)
        segments = _boundary_segments(size(vertices,2))
    end
    if isnothing(markers)
        markers = ones(UInt8,size(vertices,2))
    end
    tri.segmentlist = segments
    tri.segmentmarkerlist = markers
    if !isnothing(holes)
        tri.holelist = holes
    end
    maxarea  = Printf.@sprintf "%0.15f" h^2/2
    minangle = Printf.@sprintf "%0.15f" 30.
    (tri,_) = triangulate("pea$(maxarea)q$(minangle)Q",tri)
    F = eltype(vertices)
    P = UInt8
    I = eltype(tri.edgelist)
    HPMesh{F,I,P}(tri)
end

_boundary_segments(n) = reduce(hcat,[i,mod1(i+1,n)] for i in 1:n)

"""
    $(SIGNATURES)

builds a `HPMesh` with elements of size `h` of the rectangle `[a,b]⨱[c,d]`. If `b` and `c` are ommitted, they are assumed to be zero. 
"""
function rectmesh(a,b,c,d,h)
    vertices = Matrix{Float64}([a c;b c;b d;a d]')
    meshhp(vertices,h)
end
rectmesh(a,c,h) = rectmesh(0.,a,0.,c,h)

"""
    $(SIGNATURES)

builds a `HPMesh` with elements of size `h` of the square `[0,a]⨱[0,a]`."""
squaremesh(a,h) = rectmesh(0.,a,0.,a,h)


"""
    $(SIGNATURES)

builds a `HPMesh` with elements of size `h` of the circle with center `center` and radius `rad`. If the center is ommited, it is assumed to be the origin. If the radius is also ommited, it is assumed to be 1. 
"""
function circmesh(center,rad,h)
    cx,cy = center
    n     = Int(1+2π*rad÷h)
    θ     = range(start=0.,stop=2π,length=n)[1:end-1]
    verts = [cx .+ cos.(θ');cy .+ sin.(θ')]
    meshhp(verts,h)
end
circmesh(rad,h) = circmesh((0.,0.),rad,h)
circmesh(h) = circmesh(1.,h)


########
function Base.sizehint!(d::Dictionary,n::Int)
    sizehint!(d.indices.slots,(1+n÷8)*8)
    sizehint!(d.indices.hashes,n)
    sizehint!(d.indices.values,n)
    sizehint!(d.values,n)
end
