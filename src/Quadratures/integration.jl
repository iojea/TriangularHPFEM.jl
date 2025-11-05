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
    coefftype(integrand.op),order(integrand.op),integrand.op,m
end

function transform_matrix!(A, vert)
    @views A[:, 1] .= 0.5(vert[:, 3] - vert[:, 2])
    @views A[:, 2] .= 0.5(vert[:, 1] - vert[:, 2])
end

function transform_term!(b, vert)
    @views b .= 0.5(vert[:, 1] + vert[:, 3])
end


function integrate(::Type{Spaces.Constant},::Type{Spaces.Order{B}},op,m::Measure{M}) where {B,F,I,P,M<:HPMesh{F,I,P}}
    (;mesh,aux) = m
    degrees_of_freedom!(mesh)
    tensors = Dict{NTuple{3,P},Array{F,2*B}}()
    (;trilist,DOFs) = mesh
    (;by_tri) = DOFs
    for p in keys(by_tri)
        
    end
    ℓ = sum(x->length(x)^2,by_tri)
    J = Vector{Int32}(undef,ℓ)
    K = Vector{Int32}(undef,ℓ)
    V = Vector{Float64}(undef,ℓ)
    Aₜ = MMatrix{2,2}(zeros(2,2))
    iAₜ = MMatrix{2,2}(zeros(2,2))
    
    r = 1
    @inbounds for t in triangles(trilist)
        dofT = dof[t]
        p, pnod = pnodes(t, mesh)
        transform_matrix!(Aₜ, view(points, :, pnod))
        iAₜ .= inv(Aₜ)
        iAₜ .= iAₜ * iAₜ'
        z = vec(iAₜ)
        dAₜ = abs(det(Aₜ))
        v = zeros(dim, dim)
        for j = 1:dim, i = 1:j
            v[i, j] = S[i, j, :] ⋅ z
        end
        v = dAₜ * C' * Symmetric(v) * C
        i = repeat(dofT, dim)
        j = repeat(dofT, inner = dim)
        J[r:r+dim^2-1] = i
        K[r:r+dim^2-1] = j
        V[r:r+dim^2-1] = v
        r += dim^2
    end
    sparse(J, K, V)
end
    

################################################################################



#ref_integrate(op::OperationField) = _ref_integrate(op,FieldsType(op))
# function _ref_integrate(p::OperationField)
#     op = operation(p)
#     op in (+,-) || throw(TypeError("Only operations + and - can be integrated."))
#     q,r = arguments(p)
#     op(ref_integrate(q),ref_integrate(r))
# end