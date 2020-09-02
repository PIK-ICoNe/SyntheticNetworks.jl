export maximum_value_heuristics,cumulative_distribution_heuristics
function maximum_value_heuristics(f::Function, max_val)
    h(g,i) = f(g)[i] < max_val
    return h
end

function cumulative_distribution_heuristics(f::Function; cdf::Function)
    h(g,i) = f(g)[i] |> cdf |> x -> rand() > x
    return h
end
