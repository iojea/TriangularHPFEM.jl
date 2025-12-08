"""
    HPPoint{D,F<:Real} <: StaticArray{Tuple{D},F,1}
a `struct` for storing a coordinate point. It inherits from `StaticVector`. An `HPPoint` can be created with the coordinates or from another container:
```
   julia> HPPoint(2.0,1)
   julia> HPPoint((2.0,1))
   julia> HPPoint([2,1.0])
   julia> HPPoint(@SVector[2.0,1]) 
```
"""
struct HPPoint{D,F<:Real} <: StaticArray{Tuple{D},F,1}
    data::NTuple{D,F}
    
    function HPPoint{D,F}(x::NTuple{D,F}) where {D,F<:Real}
        #StaticArrays.check_array_parameters(Tuple{D},F,Val{1},Val{D})
        new{D,F}(x)
    end

    function HPPoint{D,F}(x::NTuple{D,Any}) where {D,F<:Real}
        #StaticArrays.check_array_parameters(Tuple{D}, F, Val{1}, Val{D})
        new{D,F}(StaticArrays.convert_ntuple(F, x))
    end
end

function HPPoint(t::Tuple)
    tt = promote(t...)
    T = eltype(tt)
    HPPoint{length(t),T}(tt...)
end

HPPoint(t::AbstractArray) = HPPoint(Tuple(t))
#HPPoint(t::StaticArray) = HPPoint(Tuple(t))
HPPoint(x,y) = HPPoint((x,y)) 
HPPoint{D,F}() where {D,F} = HPPoint(Tuple(zero(F) for _ in 1:D))

const ORIGIN = HPPoint{2,Float64}()


Base.getindex(p::HPPoint{D,F},i::I) where {D,F,I} = getfield(p,:data)[i]
Base.getindex(p::HPPoint{D,F},i::Int64) where {D,F} = getfield(p,:data)[i]


# function Base.getproperty(p::HPPoint{D,F},sym::Symbol) where {D,F}
#     if sym===:x
#         return p.data[1]
#     elseif sym===:y
#         return p.data[2]
#     elseif sym===:z
#         if D==3
#             return p.data[3]
#         else
#             throw(DomainException("The point is 2D. No `z` component is defined."))
#         end
#     else
#         return getfield(p,sym)
#     end
# end