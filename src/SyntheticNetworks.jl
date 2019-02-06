__precompile__(false)
module SyntheticNetworks

using Random
using Parameters
using LightGraphs
using SpatialIndexing

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
export SyntheticNetwork

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


end # module
