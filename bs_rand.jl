include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Statistics


function measure_b!(state::Vector{Int8}; n_layers::Integer, interaction_range::Integer, generator_set::Vector{Matrix{Rational{Int64}}}, d::Dict{Matrix{Rational{Int64}}, Vector{Vector{Int8}}})
    L = length(state)
    if interaction_range ≤ L
        state_sum = fill(Int64(0), L)
        for layer in 1:n_layers
            for sublayer in 1:interaction_range
                idx = [sublayer:sublayer+interaction_range-1;]
                while idx[end] ≤ L
                    group_elem = Int64[1 0; 0 1]
                    for g in state[idx]
                        group_elem *= generator_set[g+1]
                    end
                    state[idx] = rand(d[group_elem])
                    idx = idx .+ interaction_range
                end
            end
            state_sum = state_sum .+ (state .== 2) .- (state .== 4)
        end
    else
        raise(warning("Interaction range larger than the system size."))
    end
    return state_sum / n_layers
end


# 0 = vacuum
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

interaction_range = 3
d = deserialize("data/bs/d_L$(interaction_range).dat")

L_list = [50:10:150;]
# t_meas_list = [1, 3, 10, 30, 100, 300, 3*10^3, 10^4, 3*10^4]
# params_list = [(L, t_meas) for L in L_list, t_meas in t_meas_list]
# params = params_list[parse(Int, ARGS[1])]

# L = params[1]
# t_meas = params[2]
L = L_list[parse(Int, ARGS[1])]

# println(" L = $(L) \n t_meas = $(t_meas) \n")
println(" L = $(L) \n")

N = 10^5
# nb2 = Float64[]
nb = zeros(Float64, L)
for i in 1:N
    if mod(i, 1000) == 0
        println("i = $i")
        flush(stdout)
    end

    # Generate a random state
    state = Int8.(rand([0:length(generator_set)-1;], L))

    # In the identity sector
    while reduce(*, generator_set[state .+ 1]) ≠ e
        state = Int8.(rand([0:length(generator_set)-1;], L))
    end

    # In the 0 b-charge sector
    # while length(state[state .== 2]) ≠ length(state[state .== 4])
    #     state = Int8.(rand([0:length(generator_set)-1;], L))
    # end

    nb .+= (state .== 2) .- (state .== 4)

    # nb = measure_b!(state, n_layers=t_meas, interaction_range=interaction_range, generator_set=generator_set, d=d)
    # push!(nb2, sum(nb.^2) / L)
end
nb = nb / N

# serialize("data/bs/bs_random_words/0b_sector/bs_nb2_random_words_L$(L)_N$(int2str(N))_$(int2str(t_meas)).dat", nb2)
serialize("data/bs/bs_random_words/nb_id_sector/bs_nb_random_words_L$(L)_N$(int2str(N))_$(int2str(t_meas)).dat", nb)



