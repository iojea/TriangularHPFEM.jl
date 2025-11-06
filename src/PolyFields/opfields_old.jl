struct OperationField{F<:Function,F1,F2} <: Function
    op::F
    args::Tuple{F1,F2}
end

arguments(of::OperationField) = of.args
operation(of::OperationField) = of.op

(of::OperationField)(x,y) = operation(of)((f(x,y) for f in arguments(of))...)
(of::OperationField)(x) = of(x[1],x[2])

for op in (:+,:-,:*)
    expr = Meta.parse("(Base.:$op)(p::Function,q::Function) = OperationField($op,(p,q))")
    eval(expr)
end







