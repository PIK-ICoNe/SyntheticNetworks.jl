using SyntheticNetworks, Test
using EmbeddedGraphs: EmbeddedGraph, nv
using Random
Random.seed!(42);

const RPG = RandomPowerGrid(100, 1)
const g, t_list = generate_graph(RPG)

@testset "RPG" begin
    @test RPG isa RandomPowerGrid
    @test g isa EmbeddedGraph
    @test nv(g) == RPG.n
end

