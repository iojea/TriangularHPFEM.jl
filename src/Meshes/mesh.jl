const TriangleList{I,P,F} = Dictionary{HPTriangle{I},TriangleProperties{P,F}} where {I,P,F}
const EdgeList{I,P} = Dictionary{HPEdge{I},EdgeProperties{P,Bool}} where {I,P}
#const AuxList{P,N,F} = Dictionary{NTuple{3,P},AuxData{N,F}} where {P,N,F}

struct DOFs{I<:Integer}
    by_edge::Dictionary{HPEdge{I},Vector{I}}
    by_tri::Dictionary{HPTriangle{I},Vector{I}}
end

DOFs(::Type{I}) where I<:Integer = DOFs(Dictionary{HPEdge{I},Vector{I}}(),Dictionary{HPTriangle{I},Vector{I}}())



abstract type HPTriangulation end


"""
    $(SIGNATURES)

A mesh for `HP` finite element methods. Its fields are: 
    + `points::Vector{HPPoint{2,F}}`: a vector of points
    + `trilist::TriangleList{I,P,F}`: set of trilist
    + `edgelist::EdgeList{I,P}`: set of edgelist
    + `dofs::DOFs{I}`: auxiliary data for integrating the 

    HPMesh(tri::TriangulationIO)
builds an `HPMesh` from a `Triangulate.TriangulatioIO` struct.

    Meshes of type `HPMesh` can also be constructed using the helper functions such us `circmesh`, `rectmesh`, `squaremesh`. For a more general constructor, check the docs for `meshhp`.
"""
struct HPMesh{F<:AbstractFloat,I<:Integer,P<:Integer} <: HPTriangulation
    points::Vector{HPPoint{2,F}} #ElasticMatrix{F,Vector{F}}
    trilist::TriangleList{I,P,F}
    edgelist::EdgeList{I,P}
    dofs::DOFs{I}
end

#DEFINE THIS!
HPMesh(v,t::TriangleList{I,P,F},e::EdgeList{I,P}) where {I,P,F} = HPMesh(v,t,e,DOFs(I))

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
Given a triangle  `t`, of  a mesh `mesh` returs the degrees of the edges of `t`
"""
function degrees(t::HPTriangle{I},mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    eds = edges(t)
    p = Tuple(degree.(getindices(edgelist,eds)))
    sort(p)
end
"""
  $(SIGNATURES)

Given a triangle `t` belonging to a `mesh`, it returns the (sorted) degrees
`p₁<=p₂<=p₃` of its edges, and the edges also sorted according to the degrees.
"""
function pedges(t::HPTriangle{I},mesh::HPMesh{F,I,P}) where {F,I,P}
    p,nod = pnodes(t,mesh)
    eds   = [HPEdge(nod[SVector(i,mod1(i+1,3))]) for i in 1:3]
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

    compute_dimension(p₁,p₂,p₃)
    compute_dimension(p₁,p₂)   
    compute_dimension(p₁)
    compute_dimension(t::AbstractArray)

Computes the dimension of the space ℓp₁p₂p₃. 
"""
compute_dimension(p₁, p₂, p₃) = sum(min(p₂, p₃ - j) + 1 for j = 0:p₁);
compute_dimension(p₁, p₂) = compute_dimension(p₁, p₂, p₂)
compute_dimension(p₁) = compute_dimension(p₁, p₁)
compute_dimension(t::T) where {T<:AbstractArray} = compute_dimension(t...)

"""
    degrees_of_freedom_by_edge(mesh::HPMesh{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}

Creates a dictionary (from `Dictionaries.jl`) where the keys are the edges of `mesh` and the values are vectors with indices corresponding to the nodal degrees of freedom. 
"""
function degrees_of_freedom_by_edge!(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;points,edgelist,dofs) = mesh 
    by_edge  = dofs.by_edge
    i        = size(points,2)+1
    for edge in edges(edgelist)
        med  = collect(i:i+degree(edgelist[edge])-2)
        set!(by_edge,edge,[edge[1],med...,edge[2]])
        i   += degree(edgelist[edge])-1
    end
end
"""
    degrees_of_freedom(mesh::HPMesh{F,I,P}) where {F,I,P}

Creates a dictionary (from `Dictionaries.jl`) where the keys are the triangles of the mesh and the values are vectores storing the indices of the corresponding degrees of freedom. 

Internally, `degrees_of_freedom_by_edge` is called in order to obtain the nodal degrees of freedom, and then the degrees of freedom corresponding to bubble functions are computed. 

If a dictionary of degrees of freedom by edge has already been computed, it is recommended to run: 

    degrees_of_freedom(mesh::HPMesh{F,I,P},by_edge::Dictionary{HPEdge{I},Vector{I}}) where {F,I,P}
"""
function degrees_of_freedom!(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist,trilist,dofs) = mesh
    (;by_edge,by_tri) = dofs
    if isempty(by_edge)
        degrees_of_freedom_by_edge!(mesh)
    end
    k       = maximum(maximum.(by_edge))+1 #first non-edge dof
    for t in triangles(trilist)
        p,t_edges = pedges(t,mesh)
        newdofs = zeros(I,compute_dimension(p))
        set!(by_tri,t,newdofs)
        j = 1 #counter of dof in current triangle
        println(t_edges)
        @inbounds for i in 1:3
            println(t_edges[i])
            println(by_edge[t_edges[i]])
            newdof = by_edge[t_edges[i]]
            if  same_order(t_edges[i],edgelist)
                newdofs[j:j+length(newdof)-2] .= newdof[1:end-1]
            else
                newdofs[j:j+length(newdof)-2] .= reverse(newdof[2:end])
            end 
            j += length(newdof)-1
        end
        newdofs[j:end] = k:k+(length(newdofs)-j)
        k += length(newdofs)-j + 1
    end
end


"""
    marked_dof(mesh::HPMesh{F,I,P},marker) where {F,I,P}

Returs a list of indices corresponding to the degrees of freedom marked with `marker`. If a vector of markers is passed, it returs degrees of freedom marked with any of them. 

Internally, `marked_dof` 
"""
function marked_dof(mesh::HPMesh{F,I,P},marker::N) where {F,I,P,N<:Integer}
    marked_dof(mesh,[marker])
end
function marked_dof(mesh::HPMesh{F,I,P},marker) where {F,I,P}
    (;dofs) = mesh
    if isempty(dofs.by_edge)
        degrees_of_freedom_by_edge!(mesh)
    end
    _marked_dof(mesh,marker)
end
function _marked_dof(mesh::HPMesh{F,I,P},markerslist::AbstractVector) where {F,I,P}
    (;edgelist,dofs) = mesh
    (;by_edge) = dofs
    v = Vector{I}()
    for e in edges(edgelist)
        if marker(edgelist[e]) in markerslist
            push!(v,by_edge[e]...)
        end
    end
    return unique(v)
end

"""
    boundary_dof(mesh::HPMesh{F,I,P}) where {F,I,P}

Returns a list of the degrees of freedom lying at the boundary of the mesh.
Internally, it calls `degrees_of_freedom_by_edge` in order to obtain the nodal degrees of freedom. If a dof by edge dictionary has already been computed, it is recommended to run: 

    boundary_dof(mesh,by_edge) 
"""
function boudary_dof(mesh::HPMesh{F,I,P}) where {F,I,P}
    by_edge = degrees_of_freedom_by_edge(mesh)
    boundary_dof(mesh,by_edge)
end
function boundary_dof(mesh::HPMesh{F,I,P},by_edge) where {F,I,P}
    marked_dof(mesh,by_edge,[1,2])
end

function dirichlet_dof(mesh::HPMesh{F,I,P}) where {F,I,P}
    marked_dof(mesh,1)
end

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


# function DOFs!(mesh::HPMesh)
#     degrees_of_freedom_by_edge!(mesh)
#     degrees_of_freedom!(mesh)
# end

# """
#     degrees_of_freedom_by_edge(mesh::HPMesh{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}

# Creates a dictionary (from `Dictionaries.jl`) where the keys are the edges of `mesh` and the values are vectors with indices corresponding to the nodal degrees of freedom. 
# """
# function degrees_of_freedom_by_edge!(mesh::HPMesh{F,I,P}) where {F,I,P}
#     (;points,edgelist,dofs) = mesh 
#     by_edge  = dofs.by_edge
#     i        = size(points,2)+1
#     for edge in edges(edgelist)
#         med  = collect(i:i+degree(edgelist[edge])-2)
#         set!(by_edge,edge,SVector{degree(edgelist[edge]),I}[edge[1],med...,edge[2]])
#         i   += degree(edgelist[edge])-1
#     end
# end

# """
#     degrees_of_freedom(mesh::HPMesh{F,I,P}) where {F,I,P}

# Creates a dictionary (from `Dictionaries.jl`) where the keys are the triangles of the mesh and the values are vectores storing the indices of the corresponding degrees of freedom. 

# Internally, `degrees_of_freedom_by_edge` is called in order to obtain the nodal degrees of freedom, and then the degrees of freedom corresponding to bubble functions are computed. 

# If a dictionary of degrees of freedom by edge has already been computed, it is recommended to run: 

#     degrees_of_freedom(mesh::HPMesh{F,I,P},by_edge::Dictionary{HPEdge{I},Vector{I}}) where {F,I,P}
# """
# function degrees_of_freedom!(mesh::HPMesh{F,I,P}) where {F,I,P}
#     (;edgelist,trilist,dofs) = mesh
#     (;by_edge,by_tri) = dofs
#     if isempty(dofs.by_edge)
#         degrees_of_freedom_by_edge!(mesh)
#     end
#     k       = maximum(maximum.(by_edge))+1 #first non-edge dof
#     for t in triangles(trilist)
#         p,t_edges = pedges(t,mesh)
#         by_tri[t]  = zeros(I,compute_dimension(p))
#         j = 1 #counter of dof in current triangle
#         @inbounds for i in 1:3
#             newdof = by_edge[t_edges[i]]
#             if  same_order(t_edges[i],edgelist)
#                 by_tri[t][j:j+length(newdof)-2] .= newdof[1:end-1]
#             else
#                 by_tri[t][j:j+length(newdof)-2] .= reverse(newdof[2:end])
#             end 
#             j += length(newdof)-1
#         end
#         by_tri[t][j:end] = k:k+(length(dof[t])-j)
#         k += length(by_tri[t])-j + 1
#     end
# end


# """
#     marked_dof(mesh::HPMesh{F,I,P},marker) where {F,I,P}

# Returs a list of indices corresponding to the degrees of freedom marked with `marker`. If a vector of markers is passed, it returs degrees of freedom marked with any of them. 

# Internally, `marked_dof` 
# """
# function marked_dof(mesh::HPMesh{F,I,P},marker::N) where {F,I,P,N<:Integer}
#     marked_dof(mesh,[marker])
# end
# function marked_dof(mesh::HPMesh{F,I,P},marker) where {F,I,P}
#     (;dofs) = mesh
#     if isempty(dofs.by_edge)
#         degrees_of_freedom_by_edge!(mesh)
#     end
#     _marked_dof(mesh,marker)
# end
# function _marked_dof(mesh::HPMesh{F,I,P},markerslist::AbstractVector) where {F,I,P}
#     (;edgelist,dofs) = mesh
#     (;by_edge) = dofs
#     v = Vector{I}()
#     for e in edges(edgelist)
#         if marker(edgelist[e]) in markerslist
#             push!(v,by_edge[e]...)
#         end
#     end
#     return unique(v)
# end

# """
#     boundary_dof(mesh::HPMesh{F,I,P}) where {F,I,P}

# Returns a list of the degrees of freedom lying at the boundary of the mesh.
# Internally, it calls `degrees_of_freedom_by_edge` in order to obtain the nodal degrees of freedom. If a dof by edge dictionary has already been computed, it is recommended to run: 

#     boundary_dof(mesh,by_edge) 
# """
# function boudary_dof(mesh::HPMesh{F,I,P}) where {F,I,P}
#     by_edge = degrees_of_freedom_by_edge(mesh)
#     boundary_dof(mesh,by_edge)
# end
# function boundary_dof(mesh::HPMesh{F,I,P},by_edge) where {F,I,P}
#     marked_dof(mesh,by_edge,[1,2])
# end

# function dirichlet_dof(mesh::HPMesh{F,I,P}) where {F,I,P}
#     marked_dof(mesh,1)
# end

