include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
# using Graphs
# using Luxor
# using Karnak
# using Colors
# using Plots


# bit_lines = 3
# N_permutations = factorial(2^bit_lines)
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


# g_elems = Vector{Int8}[Int8[0:7;]]
# count = 0
# rows = Int[]
# cols = Int[]
# construct_ham = true
# while count < length(g_elems)
#     if (count % 100 == 0) || (count == length(g_elems) - 1)
#         println("g_elem: $count, found: $(length(g_elems))")
#         flush(stdout)
#     end

#     g = g_elems[count+1]
#     for gate in gate_set
#         g_new = copy(g)
#         g_new = gate.(g_new)
#         # Check if we already have such a g in our collection
#         if !(g_new in g_elems)
#             # If not, then add it to the collection
#             push!(g_elems, g_new)
#             if construct_ham
#                 # Add the corresponding matrix element
#                 push!(rows, count)
#                 push!(cols, length(g_elems) - 1)
#             end
#         elseif construct_ham
#             # If yes, then determine the index of this group element
#             ind = findfirst(x -> x == g_new, g_elems) - 1
#             # Add the corresponding matrix element
#             push!(rows, count)
#             push!(cols, ind)
#         end
#     end
#     count += 1
# end
# if construct_ham
#     # Construct the Hamiltonian matrix
#     # This sparse matrix construction will sum up coefficients with similar row and col
#     # I.e., if row[i] == row[j] and col[i] == col[j], then the matrix element in the sparse matrix will be mel[i] + mel[j]
#     ham = sparse(rows.+1, cols.+1, 1)
# end

# edges = Edge.(collect(zip(rows.+1, cols.+1)))
# g = SimpleGraph(edges)

# serialize("data/circuit/cayley_graph_S8.dat", 
#             Dict("group_elements" => g_elems,
#                  "rows" => rows,
#                  "cols" => cols,
#                  "adj_matrix" => ham,
#                  "graph" => g))

# data = deserialize("data/circuit/cayley_graph_S8.dat")
# rows = data["rows"]
# cols = data["cols"]
# adj_matrix = data["adj_matrix"]
# graph = data["graph"]

# word = [3,2,1,3,2,1]
# bits = [0:7;]
# vertex = 1
# println("permutation: $(bits), vertex: $(vertex)")
# for gate in word
#     bits = gate_set[gate+1].(bits)
#     new_vertex = findfirst(x -> x == bits, g_elems)
#     println("permutation: $(bits), vertex: $(new_vertex), edge $(vertex) -> $(new_vertex): $(has_edge(g, vertex, new_vertex))")
#     vertex = new_vertex
# end






# relations = Vector{Int8}[]
# for s in 1:6
#     push!(relations, [s,s])
# end
# for i in 1:3, j in 1:3
#     if i ≠ j
#         push!(relations, [i,j,i,j])
#     end
# end
# for i in 1:3
#     push!(relations, [i,i+3,i,i+3], [i+3,i,i+3,i])
# end
# for i in 4:6, j in 4:6
#     if i ≠ j
#         push!(relations, [i,j,i,j,i,j])
#     end
# end
# for i in 1:3, j in 4:6
#     push!(relations, [i,j,i,j,i,j,i,j], [j,i,j,i,j,i,j,i])
# end

# g_elems = Vector{Int8}[Int8[]]
# count = 0
# rows = Int[]
# cols = Int[]
# while count < length(g_elems)
#     if (count % 100 == 0) || (count == length(g_elems) - 1)
#         println("g_elem: $count, found: $(length(g_elems))")
#         flush(stdout)
#     end

#     g = g_elems[count+1]
#     for s in 1:length(gate_set)-1
#         g_new = copy(g)
#         push!(g_new, s)
#         rel_flag = false
#         for r in relations
#             if (length(r) ≤ length(g_new)) && (g_new[end-length(r)+1:end] == r)
#                 g_new = g_new[1:end-length(r)]
#                 rel_flag = true
#                 break
#             end
#         end
#         if !rel_flag
#             # If not, then add it to the collection
#             push!(g_elems, g_new)
#             # Add the corresponding matrix element
#             push!(rows, count)
#             push!(cols, length(g_elems) - 1)
#         else
#             # If yes, then determine the index of this group element
#             ind = findfirst(x -> x == g_new, g_elems) - 1
#             # Add the corresponding matrix element
#             push!(rows, count)
#             push!(cols, ind)
#         end
#     end
#     count += 1
#     if length(g_elems) > factorial(8)
#         println("Group order is too big!")
#         break
#     end
# end



# L = 6
# e_states = deserialize("data/circuit/identity_states_L$(L).dat")
# e_states = unique((states -> deleteat!(states, states .== 0)).(e_states))[2:end]
# serialize("data/circuit/identity_states_noid_L$(L).dat", e_states)




L = 12
e_states = deserialize("data/circuit/identity_states_noid_L$(L).dat")
inverses = Dict([(Int8(s), Int8(s)) for s in 0:6])
relations = Vector{Int8}[]
for state in e_states
    if !isreducible(state, relations, inverses)
        push!(relations, state)
        println(state)
        flush(stdout)
    end
end

serialize("data/circuit/relations_L$(L).dat", relations)


