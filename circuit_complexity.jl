include("src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using PyCall
using Plots
using Polynomials
using LaTeXStrings
using Statistics


# bit_lines = 3

# # Construct the gate set (a vector with functions)
# gate_set = Function[bits::Integer -> bits]    # Identity gate
# # NOT gates
# for i in 0:bit_lines-1
#     push!(gate_set, bits::Integer -> flip_bit(bits, i))
# end
# # Toffoli gates
# for i in 0:bit_lines-1, j in 0:bit_lines-1, k in j+1:bit_lines-1
#     if i ≠ j && i ≠ k
#         push!(gate_set, bits::Integer -> begin
#                                             if get_bit(bits, j) && get_bit(bits, k)
#                                                 return flip_bit(bits, i)
#                                             else
#                                                 return bits
#                                             end
#                                         end)
#     end
# end
# dof_dim = length(gate_set)


# d_comp = Dict{Vector{Int8}, Integer}()                          # complexity of each permutation sector (length of the shortest circuit realizing this permutation)
# d_comp_states = Dict{Vector{Int8}, Vector{Vector{Int8}}}()      # states that saturate the complexity (i.e., states with the number of non-identity gates equal to the complexity) for each permutation sector

# identity_perm = Int8[0:(2^bit_lines - 1);]
# d_comp[identity_perm] = 0
# d_comp_states[identity_perm] = Vector{Int8}[[]]

# Lmax = 8
# for L in 1:Lmax
#     println("L=$L -----------------------------------------------------------")
#     flush(stdout)
#     d_max_len = Dict{Vector{Int8}, Integer}()                       # number of states in each perm sector with no identity gates
#     d_max_len_states = Dict{Vector{Int8}, Vector{Vector{Int8}}}()   # states in each perm sector with no identity gates
#     for state in product(fill([1:(dof_dim-1);], L)...)
#         bits = Int8[0:2^bit_lines-1;]
#         for s in collect(state)
#             bits = map(gate_set[s+1], bits)
#         end
#         d_max_len[bits] = get(d_max_len, bits, 0) + 1
#         d_max_len_states[bits] = push!(get(d_max_len_states, bits, Vector{Int8}[]), collect(state))
#         if !(bits in keys(d_comp))
#             d_comp[bits] = L
#         end
#         if d_comp[bits] == L
#             d_comp_states[bits] = push!(get(d_comp_states, bits, Vector{Int8}[]), collect(state))
#         end
#     end
#     println("size of d_comp: ", Base.summarysize(d_comp)/1024, " Kb")
#     println("size of d_comp_states: ", Base.summarysize(d_comp_states)/(1024*1024), " Mb")
#     flush(stdout)
#     for p in keys(d_comp)
#         if !(p in keys(d_max_len))
#             d_max_len[p] = 0
#             d_max_len_states[p] = Vector{Int8}[]
#         end
#     end
#     serialize("data/circuit/d_max_len_i$L.dat", d_max_len)
#     serialize("data/circuit/d_max_len_states_i$L.dat", d_max_len_states)
# end
# serialize("data/circuit/d_comp_until_i$Lmax.dat", d_comp)
# serialize("data/circuit/d_comp_states_until_i$Lmax.dat", d_comp_states)



# # Calculate ω(K,N) - the number of N-gate circuits with complexity K.
# ω = Dict{Tuple{Int8, Int8}, Integer}()
# d_comp = deserialize("data/circuit/d_comp.dat")
# N = 12
# for i in 1:N
#     d_max_len = deserialize("data/circuit/d_max_len_i$i.dat")
#     for (p, num) in d_max_len
#         comp = d_comp[p]
#         ω[(comp, i)] = get(ω, (comp, i), 0) + num
#     end
# end
# serialize("data/circuit/omega_until_N$N.dat", ω)


# # Calculate ν(K) - the number of permutations with complexity K.
# d_comp = deserialize("data/circuit/d_comp.dat")
# ν = [count(x -> x == i, values(d_comp)) for i in 0:12]
# serialize_py("data/circuit", "nu.dat", ν)


# d_comp = deserialize("data/circuit/d_comp.dat")
# ω_list = Dict{Tuple{Int8, Int8}, Vector{<:Integer}}()
# N = 12
# for i in 1:N
#     d_max_len = deserialize("data/circuit/d_max_len_i$i.dat")
#     for (p, num) in d_max_len
#         comp = d_comp[p]
#         ω_list[(comp, i)] = push!(get(ω_list, (comp, i), Integer[]), num)
#     end
# end
# serialize_py("data/circuit", "omega_list.dat", ω_list)


N = 12
ν = deserialize("data/circuit/nu.dat")
ω_list = deserialize("data/circuit/omega_list.dat")
log_nu = log2.(ν)
S_bar_1 = Real[]
S_bar_2 = Real[]
for K in 0:N
    sorted = sort(ω_list[(K,N)])
    push!(S_bar_1, log2(mean(sorted)))
    push!(S_bar_2, log2(maximum(sorted)))
end

x_nu = [0:9;]
x_S_bar = [0:1:12;]
nu_fit = Polynomials.fit(x_nu, log_nu[x_nu .+ 1], 1)
S_bar_fit_1 = Polynomials.fit(x_S_bar, S_bar_1[x_S_bar .+ 1], 1)
S_bar_fit_2 = Polynomials.fit(x_S_bar, S_bar_2[x_S_bar .+ 1], 1)

x = [0:12;]
plot(x, log_nu, seriestype=:scatter, label=L"$\log_2(\nu)$", color=:red, windowsize=(1200, 800), legend_font=14)
plot!(x, S_bar_1, seriestype=:scatter, label=L"$\bar{S}_1$", color=:blue)
plot!(x, S_bar_2, seriestype=:scatter, label=L"$\bar{S}_2$", color=:green)
plot!(x, nu_fit.(x), label="$(nu_fit)", color=:red, lw=2)
plot!(x, S_bar_fit_1.(x), label="$(S_bar_fit_1)", color=:blue, lw=2)
plot!(x, S_bar_fit_2.(x), label="$(S_bar_fit_2)", color=:green, lw=2)
xlabel!(L"Κ")
# savefig("pics/circuit/alpha_beta")






