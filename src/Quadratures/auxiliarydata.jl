struct AuxDegData{F<:AbstractFloat,N}
    ∂nodes::Vector{SVector{2,F}}
    C::SMatrix
    local_tensor::MMatrix{N,N,F}
end

function AuxDegData{F}(p::NTuple{3,P}) where {F<:AbstractFloat,P<:Integer}
    ∂nodes = boundary_nodes(F,p)
    C = matrix_C(p,∂nodes)
    N = size(C,1)
    local_tensor = MMatrix{N,N,F}(zeros(F,N,N))
    AuxDegData{F}(∂nodes,C,local_tensor)
end

degs(::Type{P},p₁,p₂,p₃) where {P<:Integer} = (P(p₁),P(p₂),P(p₃))
degs(p₁,p₂,p₃) = degs(UInt8,p₁,p₂,p₃)

compute_dimension(p) = sum(min(p[2],p[3]-j)+1 for j in 0:p[1]);

function boundary_nodes(::Type{F},p) where {F}
    seg₁   = range(start=-one(F),stop=one(F),length=p[1]+1)
    seg₂   = range(start=-one(F),stop=one(F),length=p[2]+1)
    seg₃   = range(start=one(F),stop=-one(F),length=p[3]+1)
    L = sum(p)
    nodes  = Vector{SVector{2,F}}(undef,L)
    i = 1
    for x in seg₁[1:end-1]
        nodes[i] = @SVector[x,-one(F)]
        i += 1
    end
    for y in seg₂[1:end-1]
        nodes[i] = @SVector[one(F),y]
        i += 1
    end
    for z in seg₃[1:end-1]
        nodes[i] = @SVector[z,z]
        i += 1
    end
    nodes
end;
boundary_nodes(p) = boundary_nodes(Float64,p)

@views function matrix_F(::Type{T},p,nodes) where T
    nₙ = length(nodes)
    sb = StandardBasis(p)
    n = length(sb)
    F  = zeros(T,nₙ,n)
    for (i,b) in enumerate(sb)
        F[:,i] .= b.(nodes)
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



