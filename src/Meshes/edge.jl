struct HPEdge{I} <: HPTuple{2,I}
    data::NTuple{2,I}
end

"""
    EdgeProperties(degree::P,marker::P,refine::Bool) where P<:Integer

constructs a `struct` for storing attributes of an edge. These attributes are:
+ `degree`: degree of the polynomial approximator on the edge.
+ `marker`: a marker indicating if the edge belongs to the boundary of the domain, or to an interface or to the interior. 
+ `refine`: `true` if the edge is marked for refinement. 
"""
struct EdgeProperties{P<:Integer,Bool} 
    degree::Base.RefValue{P}
    marker::P
    refine::Base.RefValue{Bool}
    #adjacent::SVector{2,I}
end
EdgeProperties(d::P,m::P,r::Bool) where P<:Integer = EdgeProperties(Ref(d),m,Ref(r))
EdgeProperties(e::EdgeProperties)  = EdgeProperties(e.degree,e.marker,e.refine)


"""
    ismarked(e::EdgeProperties)
returs `true` if `e` is marked for refinement. 
"""
@inline ismarked(e::EdgeProperties)        = e.refine[]
"""
    degree(e::EdgeProperties)
returs the degree of `e`. 
"""
@inline degree(e::EdgeProperties)          = e.degree[]
"""
    marker(e::EdgeProperties)
returs the marker of `e`. The marker indicates if `e` is a boundary edge with Dirichlet or Neumann condition, an interior boundary, etc. 
"""
@inline marker(e::EdgeProperties)          = e.marker
"""
    mark!(e::EdgeProperties)
marks `e` for refinement.  
"""
@inline mark!(e::EdgeProperties)           = e.refine[] = true
"""
    setdegree!(e::EdgeProperties,deg)
sets the degree of `e`.  
"""
@inline setdegree!(e::EdgeProperties,deg)  = e.degree[] = deg
"""
    isinterior(e::EdgeProperties)
returns `true` if the edge is an interior one.  
"""
@inline isinterior(e::EdgeProperties)      = e.marker != 1


            
