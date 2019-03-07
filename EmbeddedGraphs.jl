
module EmbeddedGraphs

    using LightGraphs
    using Distances
    import LightGraphs.edges, Base.eltype, LightGraphs.ne, LightGraphs.nv,
            LightGraphs.has_edge, LightGraphs.has_vertex,
            LightGraphs.outneighbors, LightGraphs.vertices,
            LightGraphs.is_directed

    struct EmbeddedGraph{T1, T2} <: AbstractGraph{T1}
        graph::SimpleGraph{T1}
        vertexpos::Array{T2, 1}
        distance::Function
    end

    function Base.getindex(eg::EmbeddedGraph, i,j)
        eg.distance(eg.vertexpos[i],eg.vertexpos[j])
    end


    function distance(pos1, pos2)
        evaluate(Euclidean(), pos1, pos2)
    end

    """TODO COMMENT"""
    weights(g::EmbeddedGraph) = g

    edges(g::EmbeddedGraph, args...) = LightGraphs.edges(g.graph, args...)
    Base.eltype(g::EmbeddedGraph, args...) = Base.eltype(g.graph, args...)
    has_edge(g::EmbeddedGraph, args...) = LightGraphs.has_edge(g.graph, args...)
    has_vertex(g::EmbeddedGraph, args...) = has_vertex(g.graph, args...)
    inneighbors(g::EmbeddedGraph, args...) = inneighbors(g.graph, args...)
    ne(g::EmbeddedGraph, args...) = ne(g.graph, args...)
    nv(g::EmbeddedGraph, args...) = nv(g.graph, args...)
    outneighbors(g::EmbeddedGraph, args...) = outneighbors(g.graph, args...)
    vertices(g::EmbeddedGraph, args...) = vertices(g.graph, args...)
    is_directed(g::EmbeddedGraph, args...) = is_directed(g.graph, args...)

    function add_vertex!(g::EmbeddedGraph, pos::Array, args...)
        add_vertex!(g.graph, args...)
        push!(g.vertexpos,pos)
    end

    function rem_vertex!(g::EmbeddedGraph, vs::AbstractVector{<:Integer};
    keep_order::Bool=true, args...)
        !keep_order && @warn("keep_order needs to be true!")
        rem_vertex!(g.graph, vs, true, args...)
        deleteat!(g.vertexpos, vs)
    end

end
