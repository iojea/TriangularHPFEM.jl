struct BoundaryHPMesh{F,I,P} <:HPTriangulation
    mesh::HPMesh{F,I,P}
    marker::P
end


