struct BoundaryHPMesh{F,I,P} <:HPTriangulation
    mesh::HPMesh{F,I,P}
    marker::P
end
domainmesh(m::BoundaryHPMesh) = m.mesh
domainmesh(m::HPMesh) = m

