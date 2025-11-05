const AuxMeshData{P,F} = Dictionary{NTuple{3,P},AuxDegData{F}} where {P,F}

struct Measure{T<:HPTriangulation,P,F,Q<:QScheme}
    mesh::T
    aux::AuxMeshData{P,F}
    sch::Q
end

function Measure(mesh::HPMesh{F,P,I},degsch) where {F,P,I}
    aux = AuxMeshData{P,F}()
    for t in triangles(mesh)
        d = degrees(t,mesh)
        isin,(_,_) = gettoken(aux,d)
        if !isin
            set!(aux,d,AuxDegData(F,d))
        end
    end
    #schs = Vector{QScheme{2,F,P}}()
    #for i in round.(P,range(minim,maxim,length=n))
    #    sch = gmquadrature(Val(2),P(2i+1))
       # push!(schs,sch)
    #end
    degsch = isodd(degsch) ? P(degsch) : P(degsch+1)
    sch = gmquadrature(Val(2),degsch)
    #push!(schs,gmquadrature(Val(2),P(35)))
    Measure{typeof(mesh),P,F,QScheme{2,F,P}}(mesh,aux,sch)
end

function Measure(mesh::HPMesh{F,P,I}) where {F,P,I}
    (;edgelist) = mesh
    maxdeg = maximum(Meshes.degree.(edgelist))
    Measure(mesh,2maxdeg+1)
end



