struct Measure{T<:HPTriangulation,Q<:QScheme}
    mesh::T
    schs::Vector{Q}
end

function Measure(mesh::HPMesh{F,P,I}) where {F,P,I}
    (;trilist) = mesh
    minim,maxim = extrema(last∘degrees,trilist)
    n = 1+(maxim-minim)÷5
    schs = Vector{QScheme{2,F,P}}()
    for i in round.(P,range(minim,maxim,length=n))
        sch = gmquadrature(Val(2),P(2i+1),T̂)
        push!(schs,sch)
    end
    Measure(mesh,schs)
end


