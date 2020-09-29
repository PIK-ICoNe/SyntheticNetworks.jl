export maximum_value_heuristics,
    cumulative_distribution_heuristics, maximum_degree

function maximum_value_heuristics(f::Function, max_val)
    # h(g,i) = f(g)[i] < max_val
    return (g, i) -> f(g)[i] < max_val
end

function cumulative_distribution_heuristics(f::Function; cdf::Function)
    return (g, i) -> f(g)[i] |> cdf |> x -> rand() > x
end

function maximum_degree(k_max)
    max_degree =
        (g, i) ->
            maximum_value_heuristics(x -> degree_centrality(x, normalize = false), k_max)
    return max_degree
end

# With probability p1 + threshold, allows the connection of nodes.
# The idea is to find where V is maximal, and use this function for chosen node
# If false is returned, set V[j] = 0 and find j' where V is maximal
# Repeat this for given number of times. If none is found take the first one
# or draw new position
function accept_neighbor(neighbor_prob, type, threshold = 0)
    return rand() <= neighbor_prob[type] + threshold
end


