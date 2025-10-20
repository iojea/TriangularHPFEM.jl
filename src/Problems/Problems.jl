module Porblems
using COmmonSolve

struct HProblem{S<:AbstractSpace}
    a
    b
    space::S
    g
end
function HPProblem(a,b,space,g)
    compute_matrix(a,space)
end

function CommonSolve.solve(prob::HProblem)
    
end
end; #module
