struct DOFs{I<:Integer}
    by_edge::Dictionary{HPEdge{I},Vector{I}}
    by_tri::Dictionary{HPTriangle{I},Vector{I}}
end

DOFs(::type{I}) where I<:Integer = DOFs(Dictionary{HPEdge{I}},Vector{I}(),Dictionary{HPTriangle{I},Vector{I}}())

function DOFs!(mesh::HPMesh)
    degrees_of_freedom_by_edge!(mesh)
    degrees_of_freedom!(mesh)
end

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
        set!(by_edge,edge,SVector{degree(edgelist[edge]),I}[edge[1],med...,edge[2]])
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
    if isempty(dofs.by_edge)
        degrees_of_freedom_by_edge!(mesh)
    end
    k       = maximum(maximum.(by_edge))+1 #first non-edge dof
    for t in triangles(trilist)
        p,t_edges = pedges(t,mesh)
        by_tri[t]  = zeros(I,compute_dimension(p))
        j = 1 #counter of dof in current triangle
        @inbounds for i in 1:3
            newdof = by_edge[t_edges[i]]
            if  same_order(t_edges[i],edgelist)
                by_tri[t][j:j+length(newdof)-2] .= newdof[1:end-1]
            else
                by_tri[t][j:j+length(newdof)-2] .= reverse(newdof[2:end])
            end 
            j += length(newdof)-1
        end
        by_tri[t][j:end] = k:k+(length(dof[t])-j)
        k += length(by_tri[t])-j + 1
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

