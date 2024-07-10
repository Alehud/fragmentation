include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using Logging
using BenchmarkTools

# const ID = parse(Int, ARGS[1])

# 0 = vacuum
# 1 = a
# 2 = b
# 3 = a^(-1)
# 4 = b^(-1)

global_logger(ConsoleLogger(stderr, Logging.Error))

const interaction_range = 3
const dof_dim = 5

const L_list = [8:12;]
const L = 7

const d = deserialize("bs/data/d/d_L$(interaction_range).dat")

# function job(L::Integer)
#     # Construct Hamiltonian
#     H_terms = Vector{Tuple{Float64, Vector{Int}, Function, Function}}(undef, 0)
#     @inbounds for site in 1:(L-(interaction_range-1))
#         @views for (group_element, states) in d
#             len = length(states)
#             @inbounds for i in 1:len, j in 1:len
#                 if i ≠ j
#                     push!(H_terms, (1.0, [site:site+(interaction_range-1);], x -> x == states[i], x -> states[j]))
#                 end
#             end
#         end
#     end
#     @info "$(length(H_terms)÷2) terms in the Hamiltonian"
#     H = Hamiltonian(dof_dim, H_terms, check_hermitian=false)

#     identity_states = deserialize("bs/data/identity_states/identity_states_L$(L).dat")
#     @info "Total number of identity states: $(length(identity_states))"

#     states_all, _ = explore_full_space(H, identity_states, construct_ham=false)
#     @info "$(length(states_all)) fragile sectors found, with sizes: $(length.(states_all))"


#     @info "Serializing results...   "
#     path = "bs/data/fragile_sectors_e"
#     serialize(joinpath(mkpath(path), "fragile_sectors_e_L$(L).dat"), states_all)
#     @info "Results serialized."
# end

# @benchmark job(L)


# Construct Hamiltonian
H_terms = Vector{Tuple{Float64, Vector{Int}, Function, Function}}(undef, 0)
@inbounds for site in 1:(L-(interaction_range-1))
    @views for (group_element, states) in d
        len = length(states)
        @inbounds for i in 1:len, j in 1:len
            if i ≠ j
                push!(H_terms, (1.0, [site:site+(interaction_range-1);], x -> x == states[i], x -> states[j]))
            end
        end
    end
end
@info "$(length(H_terms)÷2) terms in the Hamiltonian"
H = Hamiltonian(dof_dim, H_terms, check_hermitian=false)

identity_states = deserialize("bs/data/identity_states/identity_states_L$(L).dat")
@info "Total number of identity states: $(length(identity_states))"


s_init = fill(Int8(0), L)
@benchmark states_all, _ = explore_connected_states(s_init, H, construct_ham=false) samples=10 seconds=100