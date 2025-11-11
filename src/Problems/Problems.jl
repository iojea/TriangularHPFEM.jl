module Porblems

using CommonSolve

struct HProblem{S<:AbstractSpace}
    a
    b
    space::S
    g
end

function HProblem(a,b,space,g)
    
end

function CommonSolve.solve(prob::HProblem)
    
end
end; #module
