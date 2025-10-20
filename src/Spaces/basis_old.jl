"""
   ```julia
    BasisIterator(::Type{P},p,xy) where P<:Integer
    ```
Creates an iterator for the Standard Basis (form my 2D Legendre Polynomials). `p` is a `Tuple` of three elements with the degrees for building the basis. Internally, the elements of `p` are converted to type `P`. If omitted, `P` is taken `UInt8`.
`xy` is a matrix of `L×2` where `xy[:,1]` and `xy[:,2]` correspond to the coordinates `x` and `y` of the points where the basis will be evaluated.
"""
struct BasisIterator{P<:Integer,T}
    p::NTuple{3,P}
    xy::Matrix{T}
    Pl::Matrix{T}
    Plm1::Matrix{T}
    buf::Vector{T}
    function BasisIterator{P,T}(p,xy) where {P,T}
        all(-1 .<= xy .<= 1) || throw(DomainError("x must lie between -1 and 1."))
        ps = issorted(p) ? p : sort(p)
        new{P,T}(maybeconvert(P,ps),maybeconvert(T,xy),ones(T,size(xy)),zeros(T,size(xy)),zeros(T,size(xy,1)))
    end
end
BasisIterator(::Type{P},p,xy) where P<:Integer = BasisIterator{P,eltype(xy)}(p,xy)
BasisIterator(p,xy) = BasisIterator(UInt8,p,xy)

maybeconvert(::Type{T},v) where T = eltype(v)==T ? v : convert.(T,v)



p₁(sb::BasisIterator) = sb.p[1]
p₂(sb::BasisIterator) = sb.p[2]
p₃(sb::BasisIterator) = sb.p[3]

Base.IteratorSize(::Type{<:BasisIterator}) = Base.HasLength()
Base.length(sb::BasisIterator) = sum(min(sb.p[2],sb.p[3]-j) for j in 0:sb.p[1]) +sb.p[1] + 1

@inline @views function Base.iterate(sb::BasisIterator{P,T}) where {P,T}
    (;Pl,Plm1) = sb
    Pl .= one(T)
    Plm1 .= zero(T)
    return copy(Pl[:,1]),(0,0)
end

function Base.iterate(sb::BasisIterator{P,T},(ℓy,ℓx)) where {P,T}
    (p₁,p₂,p₃) = sb.p
    if (ℓy == p₁) && (ℓx == min(p₃-ℓy,p₂))
        return nothing
    elseif ℓx < min(p₃-ℓy,p₂)
        return _iteratex(sb,(ℓy,ℓx))
    else
        return _iteratey(sb,(ℓy,ℓx))
    end
end

@views function _iteratex(sb::BasisIterator{P,T},(ℓy,ℓx)) where {P,T}
        (;xy,Pl,Plm1,buf) = sb
        ℓx += 1
        buf .= Pl[:,1]
        Pl[:,1] .= @. ((2ℓx-1)* xy[:,1] * buf - (ℓx-1) * Plm1[:,1])/ℓx
        Plm1[:,1] .= buf
        return Pl[:,1].*Pl[:,2],(ℓy,ℓx)
end

@views function _iteratey(sb::BasisIterator{P,T},(ℓy,ℓx)) where {P,T}
        (;xy,Pl,Plm1,buf) = sb
        ℓy += 1
        buf .= Pl[:,2]
        Pl[:,2] .= @. ((2ℓy-1)* xy[:,2] * buf - (ℓy-1) * Plm1[:,2])/ℓy
        Plm1[:,2] .= buf
        Pl[:,1] .= one(T)
        Plm1[:,1] .= zero(T)
        return copy(Pl[:,2]),(ℓy,0)
end




