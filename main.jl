include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization


# Construct U(1) toric code Hamiltonian
Lx = 2
Ly = 2
N_sites = Lx*Ly*2
H_terms = Tuple[]
for y in 0:Ly-1
    for x in 0:Lx-1
        push!(H_terms, (1.0, [(2Lx*y + 2x, '+'), (2Lx*y + 2x + 1, '+'), (2Lx*y + 2*mod(x-1, Lx), '-'), (2Lx*mod(y-1, Ly) + 2x + 1, '-')]))
        push!(H_terms, (1.0, [(2Lx*y + 2x, '+'), (2Lx*y + 2x + 1, '-'), (2Lx*y + 2*mod(x-1, Lx), '+'), (2Lx*mod(y-1, Ly) + 2x + 1, '-')]))
        push!(H_terms, (1.0, [(2Lx*y + 2x, '+'), (2Lx*y + 2x + 1, '-'), (2Lx*y + 2*mod(x-1, Lx), '-'), (2Lx*mod(y-1, Ly) + 2x + 1, '+')]))
        push!(H_terms, (1.0, [(2Lx*y + 2x, '-'), (2Lx*y + 2x + 1, '+'), (2Lx*y + 2*mod(x-1, Lx), '+'), (2Lx*mod(y-1, Ly) + 2x + 1, '-')]))
        push!(H_terms, (1.0, [(2Lx*y + 2x, '-'), (2Lx*y + 2x + 1, '+'), (2Lx*y + 2*mod(x-1, Lx), '-'), (2Lx*mod(y-1, Ly) + 2x + 1, '+')]))
        push!(H_terms, (1.0, [(2Lx*y + 2x, '-'), (2Lx*y + 2x + 1, '-'), (2Lx*y + 2*mod(x-1, Lx), '+'), (2Lx*mod(y-1, Ly) + 2x + 1, '+')]))
    end
end


# for term in H_terms
#     coef, ops = term
#     println(coef, " ", ops)
# end

s_init_arr = repeat(Bool[0, 1], Lx*Ly)
s_init = arr2int(s_init_arr)
print("s_init:  ")
print(s_init, " = ")
brint(s_init, N_sites)
states, ham, rows, cols, mels = explore_connected_states(s_init, H_terms, N_sites, construct_ham=false, check_nonzero=false, check_hermitian=false)
println("$(length(states)) states found", )
# brint(states, N_sites)
Mz = (N_sites - sum(s_init_arr))*1 + sum(s_init_arr)*(-1)
print("Mz: ", Mz)
serialize("data/Lx=$(Lx)_Ly=$(Ly)_Mz=$(Mz).dat", (length(states), states))
# print(ham)
# display(ham)
# println(rows)
# println(cols)
# println(mels)

# (N_states, states) = deserialize("Lx=$(Lx)_Ly=$(Ly)_Mz=0.dat")
# println(N_states)
# brint(states, N_sites)



