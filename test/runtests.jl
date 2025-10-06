using TriangularhpFEM
using Test

@testset "TriangularhpFEM.jl" begin
    @testset "PolyFields" include("PolyFieldsTests/runtests.jl")
end
