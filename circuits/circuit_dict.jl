include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization


bit_lines = 4
N_permutations = factorial(2^bit_lines)


# Construct the gate set (a vector with functions)
gate_set = Function[bits::Integer -> bits]    # Identity gate
# # NOT gates
# for i in 0:bit_lines-1
#     push!(gate_set, bits::Integer -> flip_bit(bits, i))
# end
# # Toffoli gates of one polarity
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
# Toffoli gates of different polarities
for polarity1 in [false, true], polarity2 in [false, true]
    for i in 0:bit_lines-1, j in 0:bit_lines-1, k in j+1:bit_lines-1
        if i ≠ j && i ≠ k
            push!(gate_set, bits::Integer -> begin
                                                if (polarity1 ⊻ get_bit(bits, j)) && (polarity2 ⊻ get_bit(bits, k))
                                                    return flip_bit(bits, i)
                                                else
                                                    return bits
                                                end
                                            end)
        end
    end
end

gateset_label = "NotTof"
foldername = "data/circuit_$(bit_lines)bits_$(gateset_label)"

for L in 1:5
    println("L: $L")
    # Construct dictionaries through exact iteration through states
    d_len, d, identity_states = construct_dict(L, gate_set, (s, gen_set) -> calc_group_elem_permutation(s, gen_set, bit_lines))
    serialize("$(foldername)/d/d_L$(L).dat", d)
    serialize("$(foldername)/d_len/d_len_L$(L).dat", d_len)
    serialize("$(foldername)/identity_states/identity_states_L$(L).dat", identity_states)

    # Construct the same dictionaries, but without any identity letters in the words
    d_noid_len, d_noid, identity_states_noid = construct_dict(L, gate_set, (s, gen_set) -> calc_group_elem_permutation(s, gen_set, bit_lines); 
                                            omit_identity_letters=true, no_identity_states=true)
    serialize("$(foldername)/d_noid/d_noid_L$L.dat", d_noid)
    serialize("$(foldername)/d_noid_len/d_noid_len_L$L.dat", d_noid_len)
    serialize("$(foldername)/identity_states_noid/identity_states_noid_L$L.dat", identity_states_noid)
    
end



# Compute d_len and d_noid_len by appending one site at a time
L_max = 100
all_perms = collect(keys(deserialize("$(foldername)/d_len/d_len_L12.dat")))
d_len = construct_d_len_iteratively(L_max, deserialize("$(foldername)/d_len/d_len_L1.dat"), (g1, g2) -> g1[g2.+1], all_perms)
serialize("$(foldername)/d_len/d_len_until_L$(L_max).dat", d_len)
d_noid_len = construct_d_len_iteratively(L_max, deserialize("$(foldername)/d_noid_len/d_noid_len_L1.dat"), (g1, g2) -> g1[g2.+1], all_perms)
serialize("$(foldername)/d_noid_len/d_noid_len_until_L$(L_max).dat", d_noid_len)


# Construct d_len by appending one site at a time
construct_d_len_iteratively((6, 10), deserialize("$(foldername)/d_len/d_len_L1.dat"), (g1, g2) -> g1[g2.+1], "$(foldername)/d_len/d_len")
# Construct d_len by appending one site at a time
construct_d_len_iteratively((6, 12), deserialize("$(foldername)/d_noid_len/d_noid_len_L1.dat"), (g1, g2) -> g1[g2.+1], "$(foldername)/d_noid_len/d_noid_len")


