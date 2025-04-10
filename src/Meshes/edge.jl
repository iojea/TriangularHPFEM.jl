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

@inline ismarked(e::EdgeProperties)        = e.refine[]
@inline degree(e::EdgeProperties)          = e.degree[]
@inline marker(e::EdgeProperties)          = e.marker
@inline mark!(e::EdgeProperties)           = e.refine[] = true
@inline setdegree!(e::EdgeProperties,deg)  = e.degree[] = deg
@inline isinterior(e::EdgeProperties)      = e.marker != 1


            
