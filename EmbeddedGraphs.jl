module EmbeddedGraphs

  using LightGraphs
  using Distances
  using SparseArrays


  import LightGraphs.edges, Base.getindex, Base.eltype, LightGraphs.ne, LightGraphs.nv,
          LightGraphs.has_edge, LightGraphs.has_vertex,
          LightGraphs.outneighbors, LightGraphs.vertices, LightGraphs.inneighbors,
          LightGraphs.is_directed, LightGraphs.edgetype, LightGraphs.weights,
          LightGraphs.add_vertex!, LightGraphs.rem_vertex!, LightGraphs.add_edge!, LightGraphs.rem_edge!

  export distance, weights, edges, EmbeddedGraph,
          has_edge, has_vertex, inneighbors, ne, nv, outneighbors, inneighbors
          vertices, is_directed, add_vertex!, rem_vertex!, edgetype,
          add_edge!, rem_edge!#, vertices_loc_x, vertices_loc_y

  struct EmbeddedGraph{T<:Integer, U} <: AbstractGraph{T}
      graph::SimpleGraph{T}
      vertexpos::Array{U, 1}
      distance::Function
  end


  function Base.getindex(eg::EmbeddedGraph, i,j)
      eg.distance(eg.vertexpos[i], eg.vertexpos[j])
  end

  function vertices_loc_x(graph::EmbeddedGraph)
      map(i->graph.vertexpos[i][1], 1:nv(graph))
  end

  function vertices_loc_y(graph::EmbeddedGraph)
      map(i->graph.vertexpos[i][2], 1:nv(graph))
  end

  function distance(pos1, pos2)
      evaluate(Euclidean(), pos1, pos2)
  end

  """TODO COMMENT"""
  function weights(g::EmbeddedGraph)
      A = spzeros(nv(g),nv(g))
      for ed in edges(g.graph)
          A[src(ed), dst(ed)] = A[dst(ed), src(ed)] = distance(g.vertexpos[src(ed)], g.vertexpos[dst(ed)])
      end
      A
  end




  edges(g::EmbeddedGraph, args...) = edges(g.graph, args...)
  Base.eltype(g::EmbeddedGraph, args...) = Base.eltype(g.graph, args...)
  has_edge(g::EmbeddedGraph, args...) = has_edge(g.graph, args...)
  has_vertex(g::EmbeddedGraph, args...) = has_vertex(g.graph, args...)
  inneighbors(g::EmbeddedGraph, args...) = inneighbors(g.graph, args...)
  ne(g::EmbeddedGraph, args...) = ne(g.graph, args...)
  nv(g::EmbeddedGraph, args...) = nv(g.graph, args...)
  outneighbors(g::EmbeddedGraph, args...) = outneighbors(g.graph, args...)
  vertices(g::EmbeddedGraph, args...) = vertices(g.graph, args...)
  is_directed(g::EmbeddedGraph, args...) = false
  is_directed(::Type{EmbeddedGraph}) = false

  is_directed(::Type{EmbeddedGraph{T, U}}) where T <: Integer where U = false
  zero(g::EmbeddedGraph) = EmbeddedGraph(SimpleGraph(0), [], distance)
  edgetype(::EmbeddedGraph{T,U}) where T <: Integer where U = LightGraphs.SimpleEdge{T}
  edgetype(g::EmbeddedGraph) = LightGraphs.SimpleEdge
  add_edge!(g::EmbeddedGraph, args...) = add_edge!(g.graph,args...)
  rem_edge!(g::EmbeddedGraph, args...) = rem_edge!(g.graph,args...)

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
