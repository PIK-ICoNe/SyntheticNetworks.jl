__precompile__(false)
module SyntheticNetworks

using Random
using LightGraphs
using SpatialIndexing

export greet

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

end # module
