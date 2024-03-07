include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization

# Symmetric group S_3 with presentation <a,b,c | a^2 = b^2 = c^3 = abc = e>
# On a square lattice (one dof on every link)
# D(w)
# D(L) = max(D(w))

# 0 = e
# 1 = a
# 2 = b
# 3 = c
# 4 = c^(-1)

e = Int64[1 0; 0 1]
a = Int64[-1 1; 0 1]
b = Int64[1 0; 1 -1]
c = Int64[-1 1; -1 0]
c̄ = Int64[0 -1; 1 -1]
generator_set = Matrix{Int64}[e, a, b, c, c̄]
generator_set_inv = (m -> Int64.(m)).(inv.(generator_set))
dof_dim = length(generator_set)

# L = 8
# # Construct dictionaries through exact iteration through states
# d_len, d, identity_states = construct_dict(L, generator_set, calc_group_elem_representation)
# serialize("data/s3/d_L$(L).dat", d)
# serialize("data/s3/d_len_L$(L).dat", d_len)
# serialize("data/s3/identity_states_L$(L).dat", identity_states)
# d_noid_len, d_noid, _ = construct_dict(L, generator_set, calc_group_elem_representation; 
#                                         omit_identity_letters=true, no_identity_states=true)
# serialize("data/s3/d_noid_L$(L).dat", d_noid)
# serialize("data/s3/d_noid_len_L$(L).dat", d_noid_len)

# d_star = Dict{NTuple{4, typeof(generator_set[1])}, Vector{Vector{Int8}}}()
# for state in product(fill([1:dof_dim;], 4)...)
#     g1 = generator_set_inv[state[2]] * generator_set_inv[state[1]]
#     g2 = generator_set_inv[state[3]] * generator_set[state[2]]
#     g3 = generator_set[state[4]] * generator_set[state[3]] 
#     g4 = generator_set[state[1]] * generator_set_inv[state[4]]
#     key = (g1, g2, g3, g4)
#     d_star[key] = push!(get(d_star, key, Vector{Int8}[]), collect(state) .- 1)
# end
# serialize("data/s3/d_star.dat", d_star)


Lx = 2
Ly = 2
d_lattice = Dict{NTuple{2, typeof(generator_set[1])}, Vector{Vector{Int8}}}()
for state in product(fill([1:dof_dim;], 2Lx*Ly)...)
    source_free = true
    for x in 1:Lx, y in 1:Ly
        if reduce(*, generator_set[collect(state[plaquette_idx(x, y; Lx=Lx, Ly=Ly)])]) ≠ e
            source_free = false
            break
        end
    end
    if source_free
        gx = reduce(*, generator_set[collect(state[row_idx(1; Lx=Lx)])])
        gy = reduce(*, generator_set[collect(state[column_idx(1; Lx=Lx, Ly=Ly)])])
        key = (gx, gy)
        d_lattice[key] = push!(get(d_lattice, key, Vector{Int8}[]), collect(state) .- 1)
    end
end


state = d_lattice[(e,e)][2]


draw_grid(state; Lx=Lx, Ly=Ly, legend=nothing)





