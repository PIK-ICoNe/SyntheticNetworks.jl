"""For the incorporation of different node types to SyntheticNetworks package."""

struct NodeType
    name::String
    probability::Float32
    neighbor_probability::Dict # Dict(Symbol(NodeType_i.name) => prob)
    method_parameter
end

NodeType() = NodeType("Node", 1., Dict(:Node => 1.), 0)

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
