"""
    QScheme

A numerical quadrature scheme. Schemes are created by calling
`gmquadrature`, and evaluated with an integration kernel function
by calling `integrate`.
"""
struct QScheme{D,T<:Number,P<:Integer}
    dim::P
    weights::Vector{T}
    points::Vector{SVector{D,T}} #check dimensions
    degree::P
end

"""
    _jump(i)
    _jump(i,dim)
A very simple function that computes recursively the jump between the number of points needed for a quadrature of degree `i`, with respect to the quadrature of degree `i-1`. Note that `i` is not the actual degree but `i=deg÷2`.
"""
@inline _jump(i) = i
@inline _jump(i,dim) = dim == 1 ? _jump(i) : sum(1:_jump(i,dim-1))

"""
   _numberofpoints(i,dim)
Computes the number of points for a quadrature of degree `deg` where `i=deg÷2`.
"""
@inline _numberofpoints(deg,dim) =  deg == 1 ? 1 : _numberofpoints(deg-1,dim)+_jump(deg,dim)

"""
   numberofpoints(deg,dim)
Computes the number of points for a quadrature of odd degree `deg`. It works for `dim=1,2,3`.
"""
numberofpoints(deg,dim)= _numberofpoints(1+deg÷2,dim)


"""
    gmquadrature(::Val{D}, degree,Tref)

# Arguments
- `D`: dimension
- `degree`: desired polynomial degree of accuracy (must be odd)
- `Tref`: simplex where the quadrature is located. A matrix of `D⨱(D+1)` with the vertices.
"""
function gmquadrature(::Val{D},degree::P,Tref) where {D,P<:Integer}
    @assert D == size(Tref,1)
    T = eltype(Tref)
    L = numberofpoints(degree,D)
    _gmquadrature(T,Val(D),Val(L),degree,Tref)
end
function gmquadrature(d::Val{D},degree::P) where {D,P}
    Tref = @SMatrix[-1. 1 1;-1. -1. 1]
    gmquadrature(d,degree,Tref)
end


_szero(::Val{D},::Type{T}) where {D,T} = MVector{D,T}(zero(T) for _ in 1:D)
_sszero(::Val{D},::Val{L},::Type{T}) where {D,L,T} = [_szero(Val(D),T) for _ in 1:L]

function _gmquadrature(::Type{T}, ::Val{D}, ::Val{L}, degree::P,Tref) where {T,D,L,P}
    D::Int
    @assert degree ≥ 0
    @assert isodd(degree)
    N = D + 1
    order = (degree - 1) ÷ 2
    exponents = get_all_exponents(Val(N), order)
    weights = zeros(T,L)
    points = _sszero(Val(D),Val(L),T)

    j = 1
    for i in 0:order
        w = T((-1)^i) * big(degree + D - 2 * i)^degree / (big(2)^(2 * order) *
             factorial(big(i)) *
             factorial(big(degree + D - i)))
        weights[j:j+length(exponents[order-i+1])-1] .= w
        for part in exponents[order - i + 1]
            points[j]  = SVector{D,T}(Tref*(map(p -> T(2 * p + 1) / (degree + D - 2 * i), part)))
            j += 1
        end
    end
    weights /= sum(weights)
    return QScheme{D,T,P}(D, weights, points, degree)
end

"""
Get all exponent combinations of dimension `dim` and maximum degree
`max_degree`. This method is actually meant for evaluating all
polynomials with these exponents. This problem is similar to the
weak_compositions, e.g.,
<https://stackoverflow.com/a/36748940/353337>. The solution here,
however, only ever adds 1, making it better suited for a possible
extension with actual polynomial evaluation.
"""
function get_all_exponents(::Val{D}, max_degree::Int) where {D}
    # Initialization, level 0
    exponents = [SVector{D,Int}(0 for d in 1:D)]

    all_exponents = Vector{SVector{D,Int}}[]
    push!(all_exponents, exponents)
    for _ in 1:max_degree
        exponents = augment(exponents)
        push!(all_exponents, exponents)
    end

    return all_exponents::Vector{Vector{SVector{D,Int}}}
end

"""
This function takes the values and exponents of a given monomial
level, e.g., [(1,0,0), (0,1,0), (0,0,1)], and augments them by one
level, i.e., [(2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)].
The method works for all dimensions and is based on the observation
that the augmentation can happen by
    1. adding 1 to all exponents from the previous level (i.e.,
       multiplication with x[0]),
    2. adding 1 to the second exponent of all exponent tuples where
       ex[0]==0 (i.e., multiplication with x[1]),
    3. adding 1 to the third exponent of all exponent tuples where
       ex[0]==0, ex[1]=0 (i.e., multiplication with x[2]),
etc. The function call is recursive.
"""
function augment(exponents::Vector{SVector{D,Int}}) where {D}
    if length(exponents) == 0 || D == 0
        return Vector{SVector{D,Int}}()
    end

    idx_leading_zero = [k for k in 1:length(exponents) if exponents[k][1] == 0]
    exponents_with_leading_zero = [popfirst(exponents[k])
                                   for k in idx_leading_zero]

    out1 = augment(exponents_with_leading_zero)
    # increment leading exponent by 1
    out = [SVector(e[1] + 1, popfirst(e)...) for e in exponents]
    append!(out, [SVector(0, e...) for e in out1])

    return out
end



# @inbounds function hatmass(p,Tref::ReferenceElement{D,N,F},C) where {D,N,F}
#     Tri,volume = Tref
#     T = eltype(Tri)
#     dim = compute_dimension(p)
#     sch  = gmquadrature(Val(2),2p[3]+1,Tref)
#     x = first.(sch.points)
#     y = getindex.(sch.points,2)
#     w = sch.weights
#     sb = BasisIterator(p,x,y)
#     n = length(sb)
#     M   = zeros(T,n,n)
#     for (j,b) in enumerate(sb)
#         for (i,c) in Iterators.filter(<=(j)∘first,enumerate(BasisIterator(p,x,y)))
#             M[i,j] += (b.*c) ⋅ w
#         end
#     end
#     SMatrix{n,n,T}(C'*Symmetric(M)*C*volume/factorial(D))
# end


# function hatstiff(p,Tref)
#     (;Tri,volume) = Tref
#     T = eltype(Tri)
#     sch  = gmquadrature(Val(2),2p[3]+1,Tref)
#     x = hyper.(first.(sch.points),0.,0.,0.)
#     dx = hyper.(first.(sch.points),1.,1.,0.)
#     y = hyper.(getindex.(sch.points,2),0.,0.,0.)
#     dy = hyper.(getindex.(sch.points,2),1.,1.,0)
#     w = sch.weights
#     n = compute_dimension(p)
#     S   = zeros(T,n,n,4)
#     for (j,b) in enumerate(BasisIterator(p,dx,y))
#         for (i,c) in  Iterators.filter(<=(j)∘first,enumerate(BasisIterator(p,dx,y)))
#             S[i,j,1] = (@. ε₁part(b)*ε₁part(c))⋅w #dxdx
#         end
#     end
#     for (j,b) in enumerate(BasisIterator(p,x,dy))
#         for (i,c) in  Iterators.filter(<=(j)∘first,enumerate(BasisIterator(p,dx,y)))
#             S[i,j,3] = (@. ε₁part(b)*ε₁part(c))⋅w #dydx
#         end
#     end
#     for (j,b) in enumerate(BasisIterator(p,dx,y))
#         for (i,c) in  Iterators.filter(<=(j)∘first,enumerate(BasisIterator(p,x,dy)))
#             S[i,j,2] = (@. ε₁part(b)*ε₁part(c))⋅w #dxdy
#         end
#     end
#     for (j,b) in enumerate(BasisIterator(p,x,dy))
#         for (i,c) in  Iterators.filter(<=(j)∘first,enumerate(BasisIterator(p,x,dy)))
#             S[i,j,4] = (@. ε₁part(b)*ε₁part(c))⋅w #dydy
#         end
#     end
#     @inbounds for k in 1:4
#         S[:,:,k] = Symmetric(S[:,:,k])
#     end
#     SArray{Tuple{n,n,4},Float64}(S*volume/factorial(D))
# end;



