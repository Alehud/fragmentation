include("src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization

# 0 = e
# 1 = a
# 2 = b
# 3 = a^(-1)
# 4 = b^(-1)

e = Int64[1 0; 0 1]
a = Int64[1 1; 0 1]
ā = Int64[1 -1; 0 1]
b = Int64[2 0; 0 1]
b̄ = [1//2 0; 0 1]
generator_set = Matrix{Rational{Int64}}[e, a, b, ā, b̄]


L = 5
# Construct dictionaries through exact iteration through states
d_len, d, identity_states = construct_dict(L, generator_set, calc_group_elem_representation)
serialize("data/bs/d_L$(L).dat", d)
serialize("data/bs/d_len_L$(L).dat", d_len)
serialize("data/bs/identity_states_L$(L).dat", identity_states)
d_noid_len, d_noid, _ = construct_dict(L, generator_set, calc_group_elem_representation; 
                                        omit_identity_letters=true, no_identity_states=true)
serialize("data/bs/d_noid_L$(L).dat", d_noid)
serialize("data/bs/d_noid_len_L$(L).dat", d_noid_len)



# Construct d_len by appending one site at a time
construct_d_len_iteratively((12, 27), deserialize("data/bs/d_len_L1.dat"), (g1, g2) -> g1*g2, "data/bs/d_len")

# Construct d_len by appending one site at a time
construct_d_len_iteratively((12, 27), deserialize("data/bs/d_noid_len_L1.dat"), (g1, g2) -> g1*g2, "data/bs/d_noid_len")




