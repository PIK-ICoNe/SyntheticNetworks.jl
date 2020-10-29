"""For the incorporation of different node types to SyntheticNetworks package."""

struct NodeType
    name::String
    probability::Float32
    neighbor_probability::Dict # Dict(Symbol(NodeType_i.name) => prob)
    method::Function # (graph::EG, vertex::Int) -> candidate::Bool
end

NodeType() = NodeType("Node", 1., Dict(:Node => 1.), default_method)
NodeType(name::String, probability::Float32, n_prob::Dict) = NodeType(name, probability, n_prob, default_method)
NodeType(name::String, probability::Float32, method::Function; n_types) = NodeType(name, probability, fill(1., n_types), default_method)
default_method(g::EmbeddedGraph,i::Int64)::Bool = true

# Step IG0
""" From a list of probabilities of drawing a node type, determines one randomly and
    returns the index of node type"""
function draw_type(node_types::Array{NodeType,1})
    prob = map(x -> x.probability, node_types)
    t0 = copy(prob)
    pushfirst!(t0 , 0.)
    type_interval = [sum(t0[1:i+1]) for i in 1:length(t0)-1]
    t_n = findfirst(type_interval .>= rand())
    return node_types[t_n]
end

function connect_types(g::EmbeddedGraph, n_types::Array{NodeType,1}, i, n)
    nodetype_frequency = countmap(n_types)
    n_prob = map(x -> x.neighbor_probability, n_types)
    weights = n_prob[Symbol(n_types[i].name)] .*
        nodetype_frequency[Symbol(n_types[i].name)]
    candidates = sample([1:i-1;i+1:nv(g)], deleteat!(weights,i), n)
    return [i in candidates for i in 1:nv(g)]
end

function connectable_nodes(g::EmbeddedGraph, n_types::Array{NodeType,1}, i)
    n_prob = map(x -> x.neighbor_probability[Symbol(n_types[i].name)], n_types)
    return n_prob .> 0.
end

