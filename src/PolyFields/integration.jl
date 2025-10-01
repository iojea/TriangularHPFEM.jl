function ref_integrate(p::PolyScalarField{F,X,N,Y,M}) where {F,X,N,Y,M}
    (;px,py) = p
    qy = integrate(py)
    x = ImmutablePolynomial((zero(F),one(F)))
    qx = px*(qy(-x)-qy(-one(F)))
    q = integrate(qx)
    q(one(F))-q(-one(F))
end

ref_integrate(op::OperationField) = _ref_integrate(op,FieldsType(op))
function _ref_integrate(p::OperationField,::PolyType)
    op = operation(p)
    op in (+,-) || throw(TypeError("Only operations + and - can be integrated."))
    q,r = arguments(p)
    op(ref_integrate(q),ref_integrate(r))
end