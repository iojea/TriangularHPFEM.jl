struct HPTriangle{I} <: HPTuple{3,I}
    data::NTuple{3,I}
end


"""
    _eval(t::HPTriangle,k)

Returns the index stored in the triangle at index `k` mod 3. 
""" 
@inline _eval(t::HPTriangle,k) = t[mod1(k,3)]


"""
    edges(t::HPTriangle)

Return a vector of edges with type `HPEdge`, containing the edges of `t`.
"""
@inline edges(t::HPTriangle)   = [HPEdge(_eval(t,i),_eval(t,i+1)) for i in 1:3]


""" 
    longestedge(t::HPTriangle)

Returns an `HPEdge` with the longest edge of `t`.
"""
 @inline longestedge(t::HPTriangle) = HPEdge(_eval(t,1),_eval(t,2))


"""
  triangle(t::T,p::AbstractMatrix)

constructs a `HPTriangle` where the vertices are given by the columns of `p[:,t]`. Hence, the triangle is defined by the indices stored in `t`, but sorted in such a way that the first edge is the longest. 
"""
function triangle(::Type{I},t,p::AbstractMatrix) where {I}
    t = I.(t)
    maxi = argmax(sum(abs2,p[:,t[SVector(1,2,3)]] - p[:,t[SVector(2,3,1)]],dims=1))[2]
    HPTriangle(t[mod1.(maxi:maxi+2,3)])
end



# Properties

"""
    TriangleProperties{P,F}(refine,η,ηₚ) where {P<:Integer,F<:AbstractFloat}

constructs a `struct` for storing attributes of a triangle. These attributes are:
+ `refine`: 
    - `0`: not marked for refinement. 
    - `1`: marked for refinement of _green_ type.
    - `2`: marked for refinement of _blue_ type.
    - `3`: marked for refinement of _red_ type. 
+ `η`: estimate for the local error. 
+ `ηₚ`: predictor of local error, based on previos estimations.  
The types can be inferred from the data:

    TriangleProperties(refine,η,ηₚ)

If only the `refine` argument is passed, `η` and `ηₚ` are initialized as `0.`
If no arguments are passed, `refine` is initialized as `Int8(0)`.
"""
struct TriangleProperties{P<:Integer,F<:AbstractFloat} 
    refine::Base.RefValue{P}
    η::Base.RefValue{F}
    ηₚ::Base.RefValue{F}
end
TriangleProperties(val::P,η::F,ηₚ::F) where {P,F} = TriangleProperties(Ref(val),Ref(η),Ref(ηₚ))
TriangleProperties(::Type{P},::Type{F}) where {P,F} = TriangleProperties(zero(P),zero(F),zero(F))

@inline ismarked(t::TriangleProperties{P,F})  where {P,F} = t.refine[] > 0
@inline isgreen(t::TriangleProperties{P,F})   where {P,F} = t.refine[] == 1
@inline isblue(t::TriangleProperties{P,F})    where {P,F} = t.refine[] ==2
@inline isred(t::TriangleProperties{P,F})     where {P,F} = t.refine[] == 3
@inline mark!(t::TriangleProperties{P,F},k)   where {P,F} = t.refine[] = k
@inline mark!(t::TriangleProperties{P,F})     where {P,F} = mark!(t,3)
@inline setη!(t::TriangleProperties{P,F},η)   where {P,F} = t.η[] = η
@inline setηₚ!(t::TriangleProperties{P,F},ηₚ) where {P,F} = t.ηₚ[] = ηₚ

