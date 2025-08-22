struct AuxDegData{N,F<:AbstractFloat}
    ∂nodes::Vector{SVector{2,F}}
    C::SMatrix{N,N,F}
    local_tensors::Vector{Array{F}}
end

function AuxDegData(::Type{F},p::NTuple{3,P}) where {F<:AbstractFloat,P<:Integer}
    ∂nodes = boundary_nodes(F,p)
    C = matrix_C(p,∂nodes)
    N = size(C,1)
    local_tensors = Vector{Array{F}}()
    AuxDegData{N,F}(∂nodes,C,local_tensors)
end

degs(::Type{P},p₁,p₂,p₃) where {P<:Integer} = (P(p₁),P(p₂),P(p₃))
degs(p₁,p₂,p₃) = degs(UInt8,p₁,p₂,p₃)

compute_dimension(p) = sum(min(p[2],p[3]-j)+1 for j in 0:p[1]);

function boundary_nodes(::Type{F},p) where {F}
    seg₁   = range(start=F(1),stop=F(-1),length=p[1]+1)
    seg₂   = range(start=F(-1),stop=F(1),length=p[2]+1)
    seg₃   = range(start=F(-1),stop=F(1),length=p[3]+1)
    L = sum(p)
    nodes  = Vector{SVector{2,F}}(undef,L)
    i = 1
    for y in seg₁[1:end-1]
        nodes[i] = @SVector[F(-1),y]
        i += 1
    end
    for x in seg₂[1:end-1]
        nodes[i] = @SVector[x,F(-1)]
        i += 1
    end
    for z in seg₃[1:end-1]
        nodes[i] = @SVector[-z,z]
        i += 1
    end
    nodes
end;
boundary_nodes(p) = boundary_nodes(Float64,p)

@views function matrix_F(::Type{T},p,nodes) where T
    xy = reinterpret(reshape,T,nodes)'
    sb = BasisIterator(p,xy)
    nₙ = size(xy,1)
    n  = length(sb)
    F  = zeros(T,nₙ,n)
    for (i,b) in enumerate(sb)
        F[:,i] .= b
    end
    return F
end;

@views function matrix_C(p,nodes)
    nₙ = length(nodes)
    F  = matrix_F(eltype(nodes[1]),p,nodes)
    n = size(F,2)  
    U,Σ,V = svd!(F,full=true)
    SMatrix{n,n}(reduce(hcat,(V[:,1:nₙ]*Diagonal(1 ./Σ)*U',V[:,nₙ+1:n])))
end



