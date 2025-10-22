const AuxMeshData{P,F} = Dictionary{NTuple{3,P},AuxDegData{F}} where {P,F}

struct Measure{T<:HPTriangulation,P,Q<:QScheme}
    mesh::T
    aux::AuxMeshData{P}
    schs::Vector{Q}
end

function Measure(mesh::HPMesh{F,P,I}) where {F,P,I}
    (;edgelist) = mesh
    minim,maxim = extrema(Meshes.degree.(edgelist))
    n = 1+(maxim-minim)÷5
    aux = AuxMeshData{P,F}()
    for t in triangles(mesh)
        d = degrees(t,mesh)
        isin,(_,_) = gettoken(aux,d)
        if !isin
            set!(aux,d,AuxDegData{F}(d))
        end
    end
    schs = Vector{QScheme{2,F,P}}()
    for i in round.(P,range(minim,maxim,length=n))
        sch = gmquadrature(Val(2),P(2i+1))
        push!(schs,sch)
    end
    push!(schs,gmquadrature(Val(2),P(35)))
    Measure{typeof(mesh),P,QScheme{2,F,P}}(mesh,aux,schs)
end


