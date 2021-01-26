__precompile__(false)

module SyntheticNetworks

import StatsBase: countmap, sample
using Random
using Parameters
using LightGraphs
using EmbeddedGraphs
using Distances
using MetaGraphs
include("NodeTypes.jl")

function __init__()
    @warn "The parameter `u` is currently not implemented."
end


"""
    SyntheticNetwork
"""

"""[1] equation
        f(i,j,G) = (d_G(i,j) + 1) ^ r / (dist_spatial(x_i,x_j))
"""
abstract type SyntheticNetwork end
export SyntheticNetwork, NodeType, RandomPowerGrid, initialise, generate_graph, grow!



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
    types::Array{NodeType,1}
    heuristic::Symbol
    function RandomPowerGrid(n, n0, p, q, r, s, u)
        new(n, n0, p, q, r, s, u)
    end
end

RandomPowerGrid(n, n0) = RandomPowerGrid(n, n0, rand(5)...,[NodeType()], :standard)
RandomPowerGrid(n,n0,p,q,r,s,u) = RandomPowerGrid(n,n0,p,q,r,s,u,[NodeType()], :standard)

function generate_graph(RPG)
    @assert RPG.heuristic in [:standard, :maxdegree, :neighbor_probability, :degree_cdf, :neighbor_probability_simple]
    # (n, n0, p, q, r, s, u) = (RPG.n, RPG.n0, RPG.p, RPG.q, RPG.r, RPG.s, RPG.s)
    eg, t_list = initialise(RPG.n0, RPG.p, RPG.q, RPG.r, RPG.s, RPG.u,
                    RPG.types, RPG.heuristic)
    grow!(eg, t_list, RPG.n, RPG.n0, RPG.p, RPG.q, RPG.r, RPG.s, RPG.u,
          RPG.types, RPG.heuristic)
    mg = Embedded_to_MetaGraph(eg, t_list)
    set_props!(mg, Dict(:pqrs => [RPG.p, RPG.q, RPG.r, RPG.s], :heuristic => RPG.heuristic))
    mg
    return mg
end


function Embedded_to_MetaGraph(eg::EmbeddedGraph, t_arr)
    mg = MetaGraph(eg.graph,1.)
    for (i,j) in enumerate(eg.vertexpos)
        set_props!(mg,i,Dict(:pos => j, :type => t_arr[i].name))
    end
    return mg
end


function initialise(n0::Int, p::Real, q::Real, r::Real, s::Real, u::Real,
                    node_types::Array{NodeType,1}, heuristic::Symbol;
                    vertex_density_prob::Function=rand_uniform_2D)
    # STEP I1
    """If the locations x_1...x_N are not given, draw them independently at
        random from ρ."""
    positions = [vertex_density_prob(i) for i=1:n0 ]
    types = [draw_type(node_types) for i = 1:n0]
    graph = EmbeddedGraph(SimpleGraph(n0), positions)

    # STEP I2
    """Initialize G to be a minimum spanning tree (MST) for x_1...x_N w.r.t.
        the distance function dist_spatial(x, y) (using Kruskal’s simple or
        Prim’s more efficient algorithm). """
    mst_graph = EmbeddedGraph(complete_graph(n0), positions)
    edges = prim_mst(mst_graph.graph, weights(mst_graph, dense=true))
    for edge in edges
        add_edge!(graph, edge)
    end

    # STEP I3
    """With probability q, draw a node i ∈ {1,...,N} uniformly at
        random, find that node l ∈ {1,...,N} which is not yet linked to i and
        for which f (i,l,G) is maximal, and add the link i–l to G."""
    m = Int(round(n0*(1-s)*(p+q), RoundDown))
    for dummy in 1:m
        i = rand(1:nv(graph))
        dist_spatial = map(j -> euclidean(graph.vertexpos[i],
            graph.vertexpos[j]), 1:nv(graph))
        #
        l_edge = Step_G34(graph, i, dist_spatial, r, types, heuristic)
        if l_edge == 0
            dummy -= 1
        else
            add_edge!(graph, l_edge, i)
        end
    end

    #"""In the new code the logic has changed and this step is equal to step G4.
    #    Put m = ⌊N_0⋅(1−s)⋅(p+q)⌋. For each a = 1...m, add a link to G as
    #    follows: Find that yet unlinked pair of distinct nodes i,j ∈ {1,...,N_0}
    #    for which f(i,j,G) [eqn. 1] is maximal, and add the link i–j to G."""
    # Step_I3!(graph, r, m) # OLD LOGIC

    return graph, types
end

function grow!(graph::EmbeddedGraph, types, n::Int, n0::Int, p, q, r, s, u,
               node_types::Array{NodeType}, heuristic::Symbol;
               vertex_density_prob::Function=rand_uniform_2D)
    for n_actual in n0+1:n
        # STEP G0
        t = draw_type(node_types)
        push!(types,t)
        """With probabilities 1−s and s, perform either steps G1–G4 or step
        G5, respectively."""
        if rand() >= s || ne(graph) == 0
            # STEP G1
            """If x i is not given, draw it at random from ρ."""
            pos = vertex_density_prob(n_actual)
            add_vertex!(graph, pos)

            # STEP G2
            """ Find that node j ∈ {1,...,N} for which dist_spatial(x_i,x_j) is
                minimal and add the link i–j to G."""
            dist_spatial = map(i -> euclidean(graph.vertexpos[nv(graph)], graph.vertexpos[i]), 1:nv(graph))
            dist_spatial[nv(graph)] = 100000. #Inf
            min_dist_vertex = argmin(dist_spatial)
            # if heuristic == :neighbor_probability_simple
            #     found = types[min_dist_vertex] !==  types[nv(graph)]
            #     while !found
            #         if rand() < 1/2
            #             dist_spatial[min_dist_vertex] = 100000.
            #             min_dist_vertex = argmin(dist_spatial)
            #             found = types[min_dist_vertex] !==  types[nv(graph)]
            #         else
            #             found = true
            #         end
            #     end
            # end
            add_edge!(graph, min_dist_vertex, nv(graph))


            # STEP G3
            """ With probability p, find that node l ∈ {1,...,N} ⍀ {j} for
                which f(i,l,G) is maximal, and add the link i–l to G."""
            if rand() <= p
                l_edge = Step_G34(graph, nv(graph), dist_spatial, r, types, heuristic)
                if l_edge !== 0
                    add_edge!(graph, l_edge, nv(graph))
                end

            end
            # STEP G4
            """ With probability q, draw a node i' ∈ {1,...,N} uniformly at
                random, find that node l' ∈ {1,...,N} which is not yet linked to
                i' and for which f(i',l',G) is maximal, and add the link i'–l'
                to G."""
            if rand() <= q
                if heuristic == :degree_cdf
                    found = false
                    while !found
                        i = rand(1:nv(graph))
                        d = degree(graph, i)
                        found = rand() <= probs(d+1,types[i].method_parameter)
                    end
                else
                    i = rand(1:nv(graph))
                end
                dist_spatial = map(j -> euclidean(graph.vertexpos[i], graph.vertexpos[j]), 1:nv(graph))
                l_edge = Step_G34(graph, i, dist_spatial, r, types, heuristic)
                if l_edge !== 0
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
            newpos = (graph.vertexpos[src(edge)] * splitval + graph.vertexpos[dst(edge)] * (1-splitval))
            add_vertex!(graph, newpos)
            add_edge!(graph, src(edge), nv(graph))
            add_edge!(graph, dst(edge), nv(graph))
            rem_edge!(graph, edge)
        end
    end
end


function Step_I3!(g::EmbeddedGraph, r::Real, m::Int)
    for dummy in 1:m
        spatial = weights(g, dense=true)
        A = floyd_warshall_shortest_paths(g, weights(g)).dists
        A = ((A .+ spatial) .^ r) ./ spatial
        map(i -> A[i,i] = 0, 1:size(A)[1])
        add_edge!(g,Tuple(argmax(A))...)
    end
end

function Step_G34(g::EmbeddedGraph, i::Int, dist_spatial, r, types::Array{NodeType}, heuristic::Symbol)
    if heuristic == :standard
        return Step_G34_standard(g, i, dist_spatial, r, types)
    elseif heuristic == :maxdegree
        return Step_G34_maxdegree(g, i, dist_spatial, r, types)
    elseif heuristic == :neighbor_probability
        return Step_G34_neighbours(g, i, dist_spatial, r, types)
    elseif heuristic == :degree_cdf
        return Step_G34_degree_cdf(g, i, dist_spatial, r, types)
    elseif heuristic == :neighbor_probability_simple
        return Step_G34_neighbours_simple(g, i, dist_spatial, r, types)
    end
end


function Step_G34_standard(g::EmbeddedGraph, i::Int, dist_spatial, r, types::Array{NodeType})
    V = dijkstra_shortest_paths(g, i).dists
    V = ((V .+ dist_spatial) .^ r) ./ dist_spatial
    V[i] = 0
    V[neighbors(g, i)] .= 0
    return argmax(V)
end


function Step_G34_maxdegree(g::EmbeddedGraph, i::Int, dist_spatial, r, types::Array{NodeType})
    d = degree_centrality(g, normalize = false)
    candidates = [d[i] .< types[i].method_parameter for i in 1:nv(g)]
    if true in candidates
        V = dijkstra_shortest_paths(g, i).dists
        V = ((V .+ dist_spatial) .^ r) ./ dist_spatial
        V[i] = 0
        V[neighbors(g, i)] .= 0
        return argmax(V[findall(candidates)])
    else
        return 0
    end
end

function Step_G34_neighbours(g::EmbeddedGraph, i::Int, dist_spatial, r, types::Array{NodeType})
    n_prob = map(x -> x.neighbor_probability[Symbol(types[i].name)], types)
    candidates = n_prob .> 0.
    if true in candidates
        V = dijkstra_shortest_paths(g, i).dists
        V = ((V .+ dist_spatial) .^ r) ./ dist_spatial
        V[i] = 0
        V[neighbors(g, i)] .= 0
        return argmax(V[findall(candidates)])
    else
        return 0
    end
end

function Step_G34_neighbours_simple(g::EmbeddedGraph, i::Int, dist_spatial, r, types::Array{NodeType})
    V = dijkstra_shortest_paths(g, i).dists
    V = ((V .+ dist_spatial) .^ r) ./ dist_spatial
    V[i] = 0
    V[neighbors(g, i)] .= 0
    j = argmax(V)
    same = types[i] == types[j]
    while same
        if rand() < 1/2
            V[j] = 0
            j = argmax(V)
            same = types[i] == types[j]
        else
            same = false
        end
    end
    return j
end


probs(k,p) = exp(p[1] + p[2] * k)
function Step_G34_degree_cdf(g::EmbeddedGraph, i::Int, dist_spatial, r, types::Array{NodeType})
    d = degree(g)
    candidates = [rand() <= probs(d[j],types[j].method_parameter) for j in 1:nv(g)]
    if true in candidates
        V = dijkstra_shortest_paths(g, i).dists
        V = ((V .+ dist_spatial) .^ r) ./ dist_spatial
        V[i] = 0
        V[neighbors(g, i)] .= 0
        return argmax(V[findall(candidates)])
    else
        return 0
    end
end

rand_uniform_2D(i) = [rand_uniform(i),rand_uniform(i)]
rand_uniform(i) = 2*(0.5-rand())


end # module
