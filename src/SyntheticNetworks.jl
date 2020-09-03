module SyntheticNetworks

using Random
using Parameters
using LightGraphs
# using SpatialIndexing
using EmbeddedGraphs
using Distances
# using MetaGraphs
import MetaGraphs: MetaGraph
include("./Heuristics.jl")


"""
    SyntheticNetwork
"""

"""[1] equation
        f(i,j,G) = (d_G(i,j) + 1) ^ r / (dist_spatial(x_i,x_j))
"""
abstract type SyntheticNetwork end
export SyntheticNetwork, RandomPowerGrid, initialise, generate_graph, f

"""
    RandomPowerGrid(num_layers, n, n0, p, q, r, s, u, sampling, α, β, γ, debug)
"""
#
# @with_kw mutable struct RandomPowerGrid <: SyntheticNetwork
#     # model parameters
#     num_layers::Int
#     n::AbstractArray{Int}
#     n0::AbstractArray{Int}
#     p::AbstractArray{Float32}
#     q::AbstractArray{Float32}
#     r::AbstractArray{Float32}
#     s::AbstractArray{Float32}
#     u::AbstractArray{Float32}
#     # sampling method
#     sampling::String
#     α::Float32
#     β::Float32
#     γ::Float32
#     # debug flag
#     debug::Bool
#     # node information (output)
#     lon::AbstractArray{Float32}
#     lat::AbstractArray{Float32}
#     lev::AbstractArray{Float32}
#
#     density::AbstractArray{Float32}
#     # internal counters
#     added_nodes::AbstractArray{Int}
#     added_edges::AbstractArray{Int}
#     n_offset::Int
#     levnodes::AbstractArray{Int}
#     cumnodes::AbstractArray{Int}
#     # internal data structures
#     levgraphs::AbstractArray
#     cumgraphs::AbstractArray
#     levrtree::AbstractArray
#     cumrtree::AbstractArray
# end
# export RandomPowerGrid
#


struct RandomPowerGrid
    n::Int
    n0::Int
    p::Float32
    q::Float32
    r::Float32
    s::Float32
    u::Float32
    t_name::Any # Array
    t_prob::Array{Float32}
    t_method::Array{Function}  # (graph::EG, vertex::Int) -> candidate::Bool
end

RandomPowerGrid(n, n0) = RandomPowerGrid(n, n0, rand(5)..., " ", [1.0], [x -> true])

function generate_graph(RPG)
    # (n, n0, p, q, r, s, u) = (RPG.n, RPG.n0, RPG.p, RPG.q, RPG.r, RPG.s, RPG.s)
    eg, t_list = initialise(
        RPG.n0,
        RPG.p,
        RPG.q,
        RPG.r,
        RPG.s,
        RPG.u,
        RPG.t_name,
        RPG.t_prob,
        RPG.t_method,
    )
    grow!(
        eg,
        t_list,
        RPG.n,
        RPG.n0,
        RPG.p,
        RPG.q,
        RPG.r,
        RPG.s,
        RPG.u,
        RPG.t_name,
        RPG.t_prob,
        RPG.t_method,
    )
    mg = MetaGraph(eg, t_list)
    return mg
    # return eg
end

# Step IG0
""" From a list of probabilities of drawing a node type, determines one randomly and
    returns the index of node type"""
function draw_type(type_prob::Array{Float32,1})
    t0 = copy(type_prob)
    pushfirst!(t0, 0.0)
    type_interval = [sum(t0[1:i+1]) for i = 1:length(t0)-1]
    typerand = rand()
    t_n = findfirst(type_interval .>= rand())
    return t_n
end

function MetaGraph(eg::EmbeddedGraph, default_weight = 1.0)
    mg = MetaGraph(eg.graph, default_weight)
    for (i, j) in enumerate(eg.vertexpos)
        set_prop!(mg, i, :pos, j)
    end
    return mg
end


function MetaGraph(eg::EmbeddedGraph, t_arr::Array{String}, default_weight = 1.0)
    mg = MetaGraph(eg.graph, default_weight)
    for (i, j) in enumerate(eg.vertexpos)
        set_props!(mg, i, Dict(:pos => j, :type => t_arr[i]))
    end
    return mg
end


function initialise(
    n0::Int,
    p::Real,
    q::Real,
    r::Real,
    s::Real,
    u::Real,
    name,
    prob::Array{Float32},
    methods::Array{Function};
    vertex_density_prob::Function = rand_uniform_2D,
)

    # we can skip this when only a single node is given
    if n0 == 1
        positions = [vertex_density_prob(n0), ]
        types = [draw_type(prob), ]
        graph = EmbeddedGraph(SimpleGraph(n0), positions)
    else
        # STEP I1
        """If the locations x_1...x_N are not given, draw them independently at
        random from ρ."""

        positions = [vertex_density_prob(i) for i = 1:n0]
        types = [draw_type(prob) for i = 1:n0]
        graph = EmbeddedGraph(SimpleGraph(n0), positions)

        # STEP I2
        """Initialize G to be a minimum spanning tree (MST) for x_1...x_N w.r.t.
            the distance function dist_spatial(x, y) (using Kruskal’s simple or
            Prim’s more efficient algorithm). """

        mst_graph = EmbeddedGraph(complete_graph(n0), positions)
        edges = prim_mst(mst_graph.graph, weights(mst_graph, dense = true))
        for edge in edges
            add_edge!(graph, edge)
        end

        # STEP I3
        """With probability q, draw a node i ∈ {1,...,N} uniformly at
            random, find that node l ∈ {1,...,N} which is not yet linked to i and
            for which f (i,l,G) is maximal, and add the link i–l to G."""

        m = Int(round(n0 * (1 - s) * (p + q), RoundDown))
        for dummy = 1:m
            i = rand(1:nv(graph))
            dist_spatial =
                map(j -> euclidean(graph.vertexpos[i], graph.vertexpos[j]), 1:nv(graph))
            #
            l_edge = Step_G34(graph, i, dist_spatial, r, methods[types])
            if l_edge == 0
                dummy -= 1
            else
                add_edge!(graph, l_edge, i)
            end
        end
    end

    #"""In the new code the logic has changed and this step is equal to step G4.
    #    Put m = ⌊N_0⋅(1−s)⋅(p+q)⌋. For each a = 1...m, add a link to G as
    #    follows: Find that yet unlinked pair of distinct nodes i,j ∈ {1,...,N_0}
    #    for which f(i,j,G) [eqn. 1] is maximal, and add the link i–j to G."""
    # Step_I3!(graph, r, m) # OLD LOGIC

    return graph, types
end

function grow!(
    graph::EmbeddedGraph,
    types,
    n::Int,
    n0::Int,
    p,
    q,
    r,
    s,
    u,
    name,
    prob::Array{Float32},
    methods::Array{Function};
    vertex_density_prob::Function = rand_uniform_2D,
)
    for n_actual = n0+1:n
        # STEP G0
        t = draw_type(prob)
        push!(types, t)
        """With probabilities 1−s and s, perform either steps G1–G4 or step
        G5, respectively."""
        if (rand() >= s) | (ne(graph) == 0)
            # STEP G1
            """If x i is not given, draw it at random from ρ."""
            pos = vertex_density_prob(n_actual)
            add_vertex!(graph, pos)

            # STEP G2
            """ Find that node j ∈ {1,...,N} for which dist_spatial(x_i,x_j) is
                minimal and add the link i–j to G."""
            dist_spatial = map(
                i -> euclidean(graph.vertexpos[nv(graph)], graph.vertexpos[i]),
                1:nv(graph),
            )
            dist_spatial[nv(graph)] = 100000.0 #Inf
            min_dist_vertex = argmin(dist_spatial)
            add_edge!(graph, min_dist_vertex, nv(graph))

            # STEP G3
            """ With probability p, find that node l ∈ {1,...,N} ⍀ {j} for
                which f(i,l,G) is maximal, and add the link i–l to G."""
            if rand() <= p
                l_edge = Step_G34(graph, nv(graph), dist_spatial, r, methods[types])
                if l_edge == 0
                    n_actual -= 1
                    continue
                else
                    add_edge!(graph, l_edge, nv(graph))
                end

            end
            # STEP G4
            """ With probability q, draw a node i' ∈ {1,...,N} uniformly at
                random, find that node l' ∈ {1,...,N} which is not yet linked to
                i' and for which f(i',l',G) is maximal, and add the link i'–l'
                to G."""
            if rand() <= q
                i = rand(1:nv(graph))
                dist_spatial =
                    map(j -> euclidean(graph.vertexpos[i], graph.vertexpos[j]), 1:nv(graph))
                l_edge = Step_G34(graph, i, dist_spatial, r, methods[types])
                if l_edge == 0
                    n_actual -= 1
                    continue
                else
                    add_edge!(graph, l_edge, i)
                end

            end

        else
            # STEP G5
            """ Select an existing link a–b uniformly at random, let
                x_i = (x_a + x_b)/2, remove the link a–b, and add two links i–a
                and i–b.
                New logic splits the nodes somewhere, not in the middle."""
            edge = rand(edges(graph))
            splitval = rand()
            newpos = (
                graph.vertexpos[src(edge)] * splitval +
                graph.vertexpos[dst(edge)] * (1 - splitval)
            )
            add_vertex!(graph, newpos)
            add_edge!(graph, src(edge), nv(graph))
            add_edge!(graph, dst(edge), nv(graph))
            rem_edge!(graph, edge)
        end
    end
end


function Step_I3!(g::EmbeddedGraph, r::Real, m::Int)
    for dummy = 1:m
        spatial = weights(g, dense = true)
        A = floyd_warshall_shortest_paths(g, weights(g)).dists
        A = ((A .+ spatial) .^ r) ./ spatial
        map(i -> A[i, i] = 0, 1:size(A)[1])
        add_edge!(g, Tuple(argmax(A))...)
    end
end


function Step_G34(g::EmbeddedGraph, i::Int, dist_spatial, r, method_arr::Array{Function,1})
    candidates = map(x -> method_arr(g, x), 1:nv(g))
    if sum(candidates) > 0
        V = dijkstra_shortest_paths(g, i).dists
        V = ((V .+ dist_spatial) .^ r) ./ dist_spatial
        V[i] = 0
        return argmax(V[candidates])
    else
        print("no node found")
        return 0
    end
end

rand_uniform_2D(i) = [rand_uniform(i), rand_uniform(i)]
rand_uniform(i) = 2 * (0.5 - rand())


end # module
