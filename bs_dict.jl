include("src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization

# 0 = vacuum
# 1 = a
# 2 = b
# 3 = a^(-1)
# 4 = b^(-1)

e = Int64[1 0; 0 1]
a = Int64[1 1; 0 1]
ā = Int64[1 -1; 0 1]
b = Int64[2 0; 0 1]
b̄ = [1//2 0; 0 1]

interaction_range = 3

# Construct the generator set (a vector with matrices)
generator_set = Matrix{Union{Int64, Rational{Int64}}}[e, a, b, ā, b̄]
dof_dim = length(generator_set)


# Construct a dictionary, where keys are group elements and values are vectors of words corresponding to these group elements
d = Dict{Matrix{Union{Int64, Rational{Int64}}}, Vector{Vector{Int8}}}()
d_len = Dict{Matrix{Union{Int64, Rational{Int64}}}, Integer}()
identity_states = Vector{Int8}[]
for state in product(fill([0:(dof_dim-1);], interaction_range)...)
    group_elem = copy(e)
    for s in collect(state)
        group_elem *= generator_set[s+1]
    end
    d[group_elem] = push!(get(d, group_elem, Vector{Int8}[]), collect(state))
    d_len[group_elem] = get(d_len, group_elem, 0) + 1
    if group_elem == e
        push!(identity_states, collect(state))
    end
end
# findfirst(states -> [0,0,0] in states, d)     # how to find a permutation realized by a specific state
# filter!(p -> length(p.second) > 1, d)   # Remove states that are not flippable
# flippable_states = collect(values(d))    # leave only values (states)

serialize("data/fragile/d_i$(interaction_range).dat", d)
serialize("data/fragile/d_len_i$(interaction_range).dat", d_len)
serialize("data/fragile/identity_states_i$(interaction_range).dat", identity_states)