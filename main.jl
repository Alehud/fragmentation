include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays


H_terms = Tuple[]
push!(H_terms, (1.0, [(0, '+'), (1, '-')]))
push!(H_terms, (1.0, [(0, '-'), (1, '+')]))
push!(H_terms, (1.0, [(1, '+'), (2, '-')]))
push!(H_terms, (1.0, [(1, '-'), (2, '+')]))
push!(H_terms, (1.0, [(2, '+'), (3, '-')]))
push!(H_terms, (1.0, [(2, '-'), (3, '+')]))
push!(H_terms, (-1.0, [(3, '+'), (0, '-')]))
push!(H_terms, (-1.0, [(3, '-'), (0, '+')]))

for term in H_terms
    coef, ops = term
    println(coef, " ", ops)
end

N_sites = 4
s_init = 4
brint(s_init, N_sites)
states, ham, rows, cols, mels = explore_connected_states(s_init, H_terms, N_sites, construct_ham=true, check_nonzero=true, check_hermitian=true)
println("States:")
brint(states)
print(ham)
display(ham)
# println(rows)
# println(cols)
# println(mels)

