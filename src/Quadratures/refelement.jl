struct ReferenceElement{D,N,F<:AbstractFloat}
    vert::SMatrix{D,N,F}
    volume::F
end

function ReferenceElement(vert::SMatrix{D,N,F}) where {D,N,F}
    N == D+1 || throw(Error("Elements other than triangles are not implemented yet"))
    vol = calc_tri_volume(vert)
    ReferenceElement(vert,vol)
end


@inbounds function calc_tri_volume(vertices::SMatrix{D,N,F}) where {D,N,F}
    X = SMatrix{D,D,F}(vertices[i, j + 1] - vertices[i, 1]
                       for i in 1:D, j in 1:D)
    vol = det(X)
    return vol
end


const T̂ = ReferenceElement(SMatrix{2,3,Float64}([-1 -1;1 -1;-1 1]),2.)

