module Porblems

using CommonSolve

struct HProblem{S<:AbstractSpace}
    a
    b
    space::S
    g
end

function HProblem(a,b,space,g)
    mock = a(space,space)
    
end

function CommonSolve.solve(prob::HProblem)
    
end
end; #module
