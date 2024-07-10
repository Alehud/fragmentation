include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using LinearAlgebra
# using Graphs
# using Luxor
# using Karnak
# using Colors
# using Plots


function nice_print(s::Vector{Int8})
    println(join(replace(s, 0 => '0', 1 => 'B', 2 => 'a', 3 => 'A')))
end

function nice_print(states::Vector{Vector{Int8}})
    for s in states
        nice_print(s)
    end
end

function num_states_in_charge_sector(L::Integer, Q::Integer)
    n = 0
    for s in product(fill(Int8[0,1,2,3], L)...)
        if count(s .== 2) - count(s .== 3) == Q
            n += 1
        end
    end
    return n
end


# 0 = 0
# B = 1
# a = 2
# A = 3
dof_dim = 4
pbc = false

s_init = Int8[2,0,1,0,3]
println("s_init: $(s_init)")


# Construct Hamiltonian
L = length(s_init)
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
for site in 1:L-1+Int(pbc)
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1], s -> (0 in s) && (1 in s), s -> Int8[1,1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1], s -> (0 in s) && ((1 in s) || (2 in s) || (3 in s)), s -> [s[2], s[1]]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1], s -> s == Int8[1,1], s -> Int8[1,0]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1], s -> s == Int8[1,1], s -> Int8[0,1]))
end
for site in 1:L-2+2*Int(pbc)
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[2,1,3], s -> Int8[0,1,0]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[0,1,0], s -> Int8[2,1,3]))
end
println(length(H_terms)÷2, " terms in the Hamiltonian")
H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);


# Explore Krylov sector
states, ham = explore_connected_states(s_init, H, construct_ham=true, verbose=false);
println("$(length(states)) states found")
# rows, cols, mels = findnz(ham)

nice_print(states)

# if isreal(ham)
#     ham = real.(ham)
# else
#     throw(error("Hamiltonian is not real."))
# end
# println("Hamiltonian: ")
# display(ham)
# ham = Matrix(ham)
# energies = eigvals(ham)
# vecs = eigvecs(ham)
# println("Eigenenergies: ", round.(energies, digits=4))
# println("Eigenvectors: ")
# for v in eachcol(vecs)
#     println(round.(v, digits=4))
# end



# for s_init in product(fill(Int8[0,1,2], L)...)
#     states, ham = explore_connected_states(collect(s_init), H, construct_ham=true, verbose=false);
#     flag = false
#     states_trim = (s -> s[s .≠ 1]).(states)
#     n = count((s -> all(s .≠ circshift(s, 1))).(states_trim))
#     if n ≠ 1
#         println(n)
#         nice_print(states)
#     end
#     # if !flag
#     #     println(s_init)
#     # end
# end