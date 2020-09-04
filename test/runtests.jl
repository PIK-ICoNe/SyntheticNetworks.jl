using SyntheticNetworks, Test
using EmbeddedGraphs: EmbeddedGraph, nv
using Random
Random.seed!(42);

const p = 0.3
const q = 0.2
const r = 1//3
const s = 0.1
const u = 0.

RPG = RandomPowerGrid(100, 1, p, q, r, s, u, " ", [1], [(g::EmbeddedGraph, i::Int) -> true])
eg, t_list = initialise(RPG.n0, RPG.p, RPG.q, RPG.r, RPG.s, RPG.u, RPG.t_name, RPG.t_prob, RPG.t_method)
grow!(eg, t_list, RPG.n, RPG.n0, RPG.p, RPG.q, RPG.r, RPG.s, RPG.u, RPG.t_name, RPG.t_prob, RPG.t_method)

@testset "RPG" begin
    @test RPG isa RandomPowerGrid
    @test eg isa EmbeddedGraph
    @test t_list isa Array{Int}
    @test nv(eg) == RPG.n
end

# you can plot an EmbeddedGraph with `using EmbeddedGraphs: gplot`.
