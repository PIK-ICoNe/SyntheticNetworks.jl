export maximum_value_heuristics,cumulative_distribution_heuristics,
    maximum_degree, maximum_betweenness
function maximum_value_heuristics(f::Function, max_val)
    # h(g,i) = f(g)[i] < max_val
    h = (g,i) -> f(g)[i] < max_val ? 1 : 0
    return h
end

function cumulative_distribution_heuristics(f::Function; cdf::Function)
    h = (g,i) -> f(g)[i] |> cdf |> x -> rand() > x ? 1 : 0
    return h
end

function maximum_degree(k_max)
    max_degree = (g,i) -> maximum_value_heuristics(
        x -> degree_centrality(x,normalize = false),
        k_max
    )
    return max_degree
end

function maximum_betweenness(b_max)
    max_betweenness = (g,i) -> maximum_value_heuristics(
        x -> betweenness_centrality(x,normalize = false),
        b_max
    )
    return max_betweenness
end

