const AuxMeshData{P,N,F} = Dictionary{NTuple{3,P},AuxDegData{N,F}} where {P,N,F}

struct Measure{T<:HPTriangulation,P,I,N,F,Q<:QScheme}
    mesh::T
    aux::AuxMeshData{P,N,F}
    schs::Vector{Q}
end

function Measure(mesh::HPMesh{F,P,I}) where {F,P,I}
    (;trilist) = mesh
    minim,maxim = extrema(last∘degrees,trilist)
    n = 1+(maxim-minim)÷5
    aux = AuxMeshData{P,N,F}()
    for t in triangles(mesh)
        d = degrees(t)
        isin,(_,_) = gettoken(mesh.trilist,e)
        if not isin
            set!(aux,d,AuxDegData(F,d))
        end
    end
    schs = Vector{QScheme{2,F,P}}()
    for i in round.(P,range(minim,maxim,length=n))
        sch = gmquadrature(Val(2),P(2i+1),T̂)
        push!(schs,sch)
    end
    Measure(mesh,aux,schs)
end


