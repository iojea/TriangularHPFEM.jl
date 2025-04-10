struct HPPoint{D,F<:AbstractFloat} <: StaticArray{Tuple{D},F,1}
    data::NTuple{D,F}
    
    function HPPoint{D,F}(x::NTuple{D,F}) where {D,F<:AbstractFloat}
        StaticArrays.check_array_parameters(Tuple{D},F,Val{1},Val{D})
        new{D,F}(x)
    end

    function HPPoint{D,F}(x::NTuple{D,Any}) where {D,F<:AbstractFloat}
        StaticArrays.check_array_parameters(Tuple{D}, F, Val{1}, Val{D})
        new{D,F}(StaticArrays.convert_ntuple(F, x))
    end
end

function HPPoint(t::Tuple)
    t = promote(t...)
    HPPoint{length(t),eltype(t)}(t)
end

HPPoint(t::AbstractArray) = HPPoint(Tuple(t))
HPPoint(t::StaticArray) = HPPoint(Tuple(t))

HPPoint{D,F}() where {D,F} = HPPoint(Tuple(zero(F) for _ in 1:D))

const ORIGIN = HPPoint{2,Float64}()


Base.getindex(p::HPPoint{D,F},i::Int) where {D,F} = getfield(p,:data)[i]

function Base.getproperty(p::HPPoint{D,F},sym::Symbol) where {D,F}
    if sym===:x
        return p.data[1]
    elseif sym===:y
        return p.data[2]
    elseif sym===:z
        if D==3
            return p.data[3]
        else
            throw(DomainException("The point is 2D. No `z` component is defined."))
        end
    else
        return getfield(p,sym)
    end
end