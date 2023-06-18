include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization


function flippable(s)
    return s[2] ≥ 2 && s[1] ≤ 1 && s[3] ≤ 1
end

function flip(s)
    return [s[1] + 1, s[2] - 2, s[3] + 1]
end

function flippable_back(s)
    return s[2] ≤ 0 && s[1] ≥ 1 && s[3] ≥ 1
end

function flip_back(s)
    return [s[1] - 1, s[2] + 2, s[3] - 1]
end


dof_dim = 3
L = 4
H_terms = Tuple[]
for i in 1:L
    push!(H_terms, (1.0, [i, mod(i+1-1, L)+1, mod(i+2-1, L)+1], flippable, flip))
    push!(H_terms, (1.0, [i, mod(i+1-1, L)+1, mod(i+2-1, L)+1], flippable_back, flip_back))
end

H = Hamiltonian(dof_dim, H_terms, check_hermitian=true)

for H_terms in H.H_terms
    println(H_terms)
end

# s_init = Int8[0, 2, 0, 2]
# println("s_init: $s_init")
# states, ham, rows, cols, mels = explore_connected_states(s_init, H, construct_ham=true)
# println("$(length(states)) states found", )
# println("states: $states")
# println("ham: ", ham)
# display(ham)
# println("rows: $rows")
# println("cols: $cols")
# println("mels: $mels")

states_all, hams = explore_full_space(H, 4, construct_ham=true)
for states in states_all
    println(states)
end




# Mz = (N_sites - sum(s_init_arr))*1 + sum(s_init_arr)*(-1)
# println("Mz: ", Mz)
# serialize("data/Lx=$(Lx)_Ly=$(Ly)_Mz=$(Mz).dat", (length(states), states))
# (N_states, states) = deserialize("Lx=$(Lx)_Ly=$(Ly)_Mz=0.dat")
# println(N_states)
# brint(states, N_sites)



