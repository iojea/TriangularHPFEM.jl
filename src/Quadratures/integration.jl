function ref_integrate(p::TensorPolynomial{F,X,N,Y,M}) where {F,X,N,Y,M}
    (;px,py) = p
    qy = Polynomials.integrate(py)
    x = ImmutablePolynomial((zero(F),one(F)))
    qx = px*(qy(x)-qy(-one(F)))
    q = Polynomials.integrate(qx)
    q(one(F))-q(-one(F))
end

ref_integrate(p::PolySum) = ref_integrate(p.left) + ref_integrate(p.right)
# function ref_integrate(p::BivariatePolynomial{F,X,N,Y,M}) where {F,X,N,Y,M}
#     x = ImmutablePolynomial((zero(F),one(F)),X)
#     ip = integrate(p.p)
#     q = ip(x) - ip(-one(F))
#     iq = integrate(q)
#     iq(one(F))-iq(-one(F))
# end

"""
    integrate(fun, scheme)
    integrate(fun, scheme, vertices::AbstractVector)
    integrate(fun, scheme, vertices::AbstractMatrix)

# Arguments
- `fun`: integrand, should accept an `SVector` as argument
- `scheme`: quadrature scheme
- `vertices`: vertices of the simplex

The vertices need to be passed either as a vector-of-vectors or as a
matrix. In the first case, there need to be `D+1` points with `D`
coordinates each. In the second case, the matrix needs to have size
`D`×`D+1`.

If the vertices are omitted, the function is called with barycentric
coordinates instead.
"""

@inbounds function ref_integrate(fun, scheme::QScheme{N,T}) where {N,T}
    @assert N > 0

    ws = scheme.weights
    ps = scheme.points
    @assert length(ws) == length(ps)

    p1 = ps[1]
    R = typeof(ws[1] * fun(p1))

    s = zero(R)
    @simd for i in 1:length(ws)
        w = ws[i]
        p = ps[i]
        s += w * fun(p)
    end

    return 2s / factorial(N - 1)
end

@inbounds function integrate(fun, scheme::QScheme{N,T},
                             vertices::SMatrix{D,N,U}) where {N,T,D,U}
    @assert N > 0
    @assert N >= D + 1

    ws = scheme.weights
    ps = scheme.points
    @assert length(ws) == length(ps)

    p1 = ps[1]
    x1 = (vertices * p1)::SVector{D}
    X = typeof(x1)
    R = typeof(ws[1] * fun(x1))

    s = zero(R)
    @simd for i in 1:length(ws)
        w = ws[i]
        p = ps[i]
        x = vertices * p
        s += w * fun(x)
    end

    # If `U` is an integer type, then Linearalgebra.det converts to
    # floating-point values; we might want a different type
    vol = R(calc_vol(vertices)) / factorial(N - 1)
    return vol * s
end
function integrate(fun, scheme::QScheme, vertices::SMatrix)
    return error("Wrong dimension for vertices matrix")
end
@inbounds function integrate(fun, scheme::QScheme{N,T},
                             vertices::SVector{N,SVector{D,U}}) where {N,T,D,U}
    return integrate(fun, scheme,
                     SMatrix{D,N,U}(vertices[n][d] for d in 1:D, n in 1:N))
end
function integrate(fun, scheme, vertices::SVector{N,<:SVector}) where {N}
    return error("Wrong dimension for vertices array")
end
function integrate(kernel, scheme::QScheme{N},
                   vertices::AbstractVector) where {N}
    @assert length(vertices) == N
    @assert N > 0
    D = length(vertices[1])
    @assert N >= D + 1
    vertices′ = SVector{N}(map(SVector{D}, vertices))
    return integrate(kernel, scheme, vertices′)
end
function integrate(kernel, scheme::QScheme{N},
                   vertices::AbstractMatrix) where {N}
    @assert size(vertices, 1) == N
    @assert N > 0
    D = size(vertices, 2)
    @assert N >= D + 1
    vertices′ = SMatrix{N,D}(vertices)'
    return integrate(kernel, scheme, vertices′)
end


function Base.:*(integrand::Integrand,m::Measure)
    integrate(coefftype(integrand.op),order(integrand.op),integrand.op,m)
end

function integrate(::Type{Spaces.Constant},::Type{Spaces.Order{B}},op,m::Measure) where B
    println("yes")
end


################################################################################



#ref_integrate(op::OperationField) = _ref_integrate(op,FieldsType(op))
# function _ref_integrate(p::OperationField)
#     op = operation(p)
#     op in (+,-) || throw(TypeError("Only operations + and - can be integrated."))
#     q,r = arguments(p)
#     op(ref_integrate(q),ref_integrate(r))
# end