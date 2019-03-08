__precompile__(false)

module SyntheticNetworks

using Random
using Parameters
using LightGraphs
# using SpatialIndexing
using EmbeddedGraphs
"""
    greet()

A function.
Here's some inline maths: ``\\sqrt[n]{1 + x + x^2 + \\ldots}``.

Here's an equation:

``\\frac{n!}{k!(n - k)!} = \\binom{n}{k}``

```jldoctest
a = 1
b = 2
a + b

# output

3
```
"""
greet() = println("Hello World!")
export greet

"""
    SyntheticNetwork
"""
abstract type SyntheticNetwork end
export SyntheticNetwork, Grid_Growth_Parameters, initialise, generate_graph, f

"""
    RandomPowerGrid(num_layers, n, n0, p, q, r, s, u, sampling, α, β, γ, debug)
"""
@with_kw mutable struct RandomPowerGrid <: SyntheticNetwork
    # model parameters
    num_layers::Int
    n::AbstractArray{Int}
    n0::AbstractArray{Int}
    p::AbstractArray{Float32}
    q::AbstractArray{Float32}
    r::AbstractArray{Float32}
    s::AbstractArray{Float32}
    u::AbstractArray{Float32}
    # sampling method
    sampling::String
    α::Float32
    β::Float32
    γ::Float32
    # debug flag
    debug::Bool
    # node information (output)
    lon::AbstractArray{Float32}
    lat::AbstractArray{Float32}
    lev::AbstractArray{Float32}

    density::AbstractArray{Float32}
    # internal counters
    added_nodes::AbstractArray{Int}
    added_edges::AbstractArray{Int}
    n_offset::Int
    levnodes::AbstractArray{Int}
    cumnodes::AbstractArray{Int}
    # internal data structures
    levgraphs::AbstractArray
    cumgraphs::AbstractArray
    levrtree::AbstractArray
    cumrtree::AbstractArray
end
export RandomPowerGrid



struct Grid_Growth_Parameters
    n::Int
    n0::Int
    p::Float32
    q::Float32
    r::Float32
    s::Float32
    u::Float32
end

Grid_Growth_Parameters(n, n0, args...) = Grid_Growth_Parameters(n0, zeros(5)...)

function generate_graph(gp)
    n = gp.n
    n0 = gp.n0
    p = gp.p
    q = gp.q
    r = gp.r
    s = gp.s
    u = gp.u
    g = initialise(n0, p, q, r, s, u)
    grow!(g, n, n0, p, q, r, s, u)
    g
end

function grow!(g::EmbeddedGraph, n::Int, n0::Int, p, q, r, s, u; vertex_density_prob=_rand_uniform)
    for dummy in n0+1:n
        if rand() >= s
            # STEP G1
            pos = [vertex_density_prob(), vertex_density_prob()]
            add_vertex!(g, pos)

            # STEP G2
            d_spatial_i = map(i -> distance(g.vertexpos[nv(g)], g.vertexpos[i]), 1:nv(g))
            d_spatial_i[nv(g)] = 100 #Inf
            min_dist_vertex = argmin(d_spatial_i)
            add_edge!(g, min_dist_vertex, nv(g))

            # STEP G3
            if rand() <= p
                l_edge = _G34(g, nv(g), d_spatial_i, r)
                add_edge!(g, l_edge, nv(g))
            end
            # STEP G4
            if rand() <= q
                i = rand(1:nv(g))
                d_spatial_i = map(j -> distance(g.vertexpos[i], g.vertexpos[j]), 1:nv(g))
                l_edge = _G34(g, nv(g), d_spatial_i, r)
                add_edge!(g, l_edge, i)
            end

        else
            # STEP G5
            edge = rand(edges(g))
            newpos = (g.vertexpos[src(edge)] + g.vertexpos[dst(edge)]) / 2.
            add_vertex!(g, newpos)
            add_edge!(g, src(edge), nv(g))
            add_edge!(g, dst(edge), nv(g))
            rem_edge!(g, edge)
        end
    end
end

function initialise(n0::Int, p, q, r, s, u; vertex_density_prob=_rand_uniform)
    # STEP I1
    positions = [[vertex_density_prob(), vertex_density_prob()] for i=1:n0 ]
    graph = EmbeddedGraph(SimpleGraph(n0), positions)
    # STEP I2
    mst_graph = EmbeddedGraph(CompleteGraph(n0), positions)
    edges = prim_mst(mst_graph, weights(mst_graph, dense=true))
    for i in edges
        add_edge!(graph,i)
    end
    # STEP I3
    m = Int(round(n0*(1-s)*(p+q), RoundDown))
    _I3!(graph, r, m)
    graph
end

function _I3!(g::EmbeddedGraph, r::Real, m::Int)
    for dummy in 1:m
        spatial = weights(g, dense=true)
        A = floyd_warshall_shortest_paths(g, weights(g)).dists
        A = ((A .+ spatial) .^ r) ./ spatial
        map(i -> A[i,i] = 0, 1:size(A)[1])
        add_edge!(g,Tuple(argmax(A))...)
    end
end

function _G34(g::EmbeddedGraph, i::Int, d_spatial_i, r)
    V = dijkstra_shortest_paths(g, i).dists
    V = ((V .+ d_spatial_i) .^ r) ./ d_spatial_i
    V[i] = 0
    argmax(V)
end

function _rand_uniform()
    2*(0.5-rand())
end

end # module
