
struct Operation{F<:Function,S1,S2}
    operator::F
    left::S1
    right::S2
end
Operation(a,c,b)  = Operation{typeof(a),typeof(b),typeof(c)}(a,b,c)
Operation(a,b) = Operation{typeof(a),typeof(b),Nothing}(a,b,nothing)

const Sp = Union{AbstractSpace,Operation}

Base.:*(n::Number,op::Sp) = Operation(*,n,op)
Base.:*(A::AbstractArray,op::Sp) = Operation(*,A,op)
Base.:*(a::Sp,b::Sp) = Operation(*,a,b)

LinearAlgebra.dot(a::AbstractArray,op::Sp) = Operation(dot,a,op)
LinearAlgebra.dot(a::Sp,b::Sp) = Operation(dot,a,b)

coefftype(op::Operation) = promote_type(coefftype(op.left),coefftype(op.right))

order(op::Operation) = order(op.left)+order(op.right)

