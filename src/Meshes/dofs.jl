struct DOFs{I<:Integer}
    by_edge::Dictionary{HPEdge{I},Vector{I}}
    by_tri::Dictionary{HPTriangle{I},Vector{I}}
end

DOFs(::Type{I}) where I<:Integer = DOFs(Dictionary{HPEdge{I}},Vector{I}(),Dictionary{HPTriangle{I},Vector{I}}())

