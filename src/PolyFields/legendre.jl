"""
```
    LegendreIterator(N)
```
Builds an iterator over de Legendre polynomials.
```
    julia> for p in LegendreIterator(4)
           println(p)
           end
    1.0
    1.0*x
    -0.5 + 1.5*x^2
    -1.5*x + 2.5*x^3
```
"""
struct LegendreIterator{I<:Integer,F<:Number,X}
    N::I
    LegendreIterator{I,F,X}(N) where {I,F,X} = new{I,F,X}(I(N))
end
LegendreIterator(N::I) where {I<:Integer} = LegendreIterator{I,Float64,:x}(N)
Base.IteratorSize(::Type{<:LegendreIterator}) = Base.HasLength()
Base.length(l::LegendreIterator{I,F,X}) where {I,F,X} = l.N

function Base.iterate(l::LegendreIterator{I,F,X}) where {I,F,X}
    q = ImmutablePolynomial((F(1)),X)
    z = zero(q)
    q,(1,q,z)
end

function Base.iterate(l::LegendreIterator{I,F,X},state) where {I,F,X}
    if state[1] == l.N
        return nothing
    else
        _iterate(l,state)
    end
end

function _iterate(::LegendreIterator{I,F,X},state) where {I,F,X}
    n,p,pm = state
    if n==1
        q = ImmutablePolynomial((zero(F),one(F)),X)
        return q,(2,q,p)
    else
        q = ((2n-1)ImmutablePolynomial((zero(F),one(F)),X)*p - (n-1)*pm)/n
        return q,(n+1,q,p)
    end
end

###########################
#     STANDARD BASIS
###########################


struct StandardBasis{P<:Integer,F<:Number,X,Y}
    degs::NTuple{3,P}
    Lx::LegendreIterator{P,F,X}
    Ly::LegendreIterator{P,F,Y}
    function StandardBasis{P,F,X,Y}(p) where {P,F,X,Y}
        p = P.(p)
        p[1]+p[2] >= p[3] || throw(ArgumentError("Degrees does not satisfy p conformity."))
        Lx = LegendreIterator{P,F,X}(p[2])
        Ly = LegendreIterator{P,F,Y}(p[1])
        new{P,F,X,Y}(p,Lx,Ly)
    end
end
StandardBasis(p::NTuple{3,P}) where P = StandardBasis{P,Float64,:x,:y}(p)
StandardBasis(p₁::P,p₂::P,p₃::P) where P = StandardBasis((p₁,p₂,p₃))
Base.IteratorSize(::Type{<:StandardBasis}) = Base.HasLength()
Base.length(sb::StandardBasis) = sum(min(sb.p[2],sb.p[3]-j) for j in 0:sb.p[1]) +sb.p[1] + 1

function Base.iterate(sb::StandardBasis{P,F,X,Y}) where {P,F,X,Y}
    (;Lx,Ly) = sb
    px,stx = iterate(Lx)
    py,sty = iterate(Ly)
    PolyScalarField(px,py),(stx,sty)
end

function Base.iterate(sb::StandardBasis{P,F,X,Y},state) where {P,F,X,Y}
    (;degs,Lx,Ly) = sb
    (p₁,p₂,p₃) = degs
    stx,sty = state
    if (sty[1] == p₁) && (stx[1] == min(p₃-sty[1],p₂))
        return nothing
    elseif stx[1] < min(p₃-sty[1],p₂)
        px,stxnew = iterate(Lx,stx)
        return PolyScalarField(px,sty[2]),(stxnew,sty)
    else
        py,stynew = iterate(Ly,sty)
        px,stxnew = iterate(Lx)
        return PolyScalarField(px,py),(stxnew,stynew)
    end
end


