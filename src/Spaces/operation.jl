
struct Operation{F<:Function,S1,S2}
    operator::F
    left::S1
    right::S2
end
Operation(a,c,b)  = Operation{typeof(a),typeof(b),typeof(c)}(a,b,c)
Operation(a,b) = Operation{typeof(a),typeof(b),Nothing}(a,b,nothing)


Base.:*(A::AbstractArray,op::Operation) = Operation(*,A,op)
Base.:*(a::AbstractSpace,b::AbstractSpace) = Operation(*,a,b)

LinearAlgebra.dot(a::AbstractSpace,b::AbstractSpace) = Operation(dot,a,b)

