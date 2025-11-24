
struct FEProblem{S<:AbstractSpace}
    a::Form
    b::Form
    space::S
    g
end

function FEProblem(a,b,space,g)
    matrix = integrate(Val(2),a,space)
    #rhs = integrate(Val(1),b,space)
end

function CommonSolve.solve(prob::FEProblem)
    nothing
end
