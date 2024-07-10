include("../src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Statistics
# using StatsBase


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

ID = parse(Int, ARGS[1])

L_list = [50:10:130;]
t_meas_list = [30, 100, 300, 10^3, 3*10^3, 10^4, 3*10^4, 10^5, 3*10^5]
experiments_list = [3*10^5, 10^5, 3*10^4, 10^4, 3*10^3, 10^3, 300, 100, 30]

L = L_list[ID]
n = L÷10
n_waves = 1
interaction_range = 3

t_meas = t_meas_list[ID]
experiments = experiments_list[ID]



println(" L = $(L) \n n = $n \n n_waves = $(n_waves) \n t_meas = $(t_meas) \n experiments = $(experiments) \n")

d = deserialize("bs/data/d/d_L$(interaction_range).dat")

# nb2_distribution = deserialize("bs/data/bs_random_words/0b_sector/bs_nb2_random_words_L$(L)_N10^5_$(int2str(t_meas)).dat")
# h = StatsBase.fit(Histogram, nb2_distribution)
# nb2_threshold = (h.edges[1][argmax(h.weights)] + h.edges[1][argmax(h.weights)+1])/2
# println(" nb2_threshold = $(nb2_threshold) \n")

tt = Int64[]
# evals_all = Vector{Vector{Float64}}[]
# nb2_all = Vector{Float64}[]

# roll_window = 100
therm_fraction = 0.1
t_after_therm = 100
print_rate = experiments÷1000 + Int(experiments < 1000)
for exper in 1:experiments
    print_flag = mod(exper, print_rate) == 0
    if print_flag
        println("----------------------------------------------------------------------------")
        println("Experiment: $exper")
        flush(stdout)
    end

    # All e state
    # state = fill(Int8(0), L)

    # w_large word sprinkled with identities
    state = repeat(Int8[fill(2, n); 1; fill(4, n); 3; fill(2, n); 3; fill(4, n); 1], n_waves)
    while length(state) < L
        insert!(state, rand(1:length(state)), 0)
    end

    # density wave of b's with other letters random
    # state = repeat(Int8[fill(2, n); fill(4, n); fill(2, n); fill(4, n)], n_waves)
    # while length(state) < L
    #     insert!(state, rand(1:length(state)), rand([0,1,3]))
    # end

    # Random state in the identity sector
    # state = Int8.(rand([0:length(generator_set)-1;], L))
    # while reduce(*, generator_set[state .+ 1]) ≠ e
    #     state = Int8.(rand([0:length(generator_set)-1;], L))
    # end

    # init_state = copy(state)
    if print_flag
        println("initial state: $state")
        flush(stdout)
    end

    # evals = Vector{Float64}[]
    nb2 = Float64[]
    
    t = 1
    t_stop = 0
    therm_flag = false
    while true
        if print_flag
            if mod(t, 100) == 0
                println("t = $t")
                flush(stdout)
            end
        end

        nb = measure_b!(state, n_layers=t_meas, interaction_range=interaction_range, generator_set=generator_set, d=d)
        # push!(evals, nb)
        push!(nb2, sum(nb.^2)/L)

        # Stop after maximal time
        if t > 5000
            push!(tt, -1)
            if print_flag
                println("Experiment did not thermalize")
                flush(stdout)
            end
            break
        end
        # Stop if the value of nb2 reaches some value
        if therm_flag == false
            if nb2[end] ≤ nb2[1] * therm_fraction
                push!(tt, t*t_meas)
                therm_flag = true
                if print_flag
                    println("Thermalized at t=$(t)")
                    flush(stdout)
                end
            end
        end
        # Stop if the rolling average of nb2 reaches some value
        # if therm_flag == false && t > roll_window
        #     if mean(nb2[end-roll_window:end]) ≤ mean(nb2[1:roll_window]) * therm_fraction
        #         push!(tt, t*t_meas)
        #         therm_flag = true
        #         println("Thermalized at t=$(t) (with roll_window=$(roll_window))")
        #         flush(stdout)
        #     end
        # end
        # If the stop condition has been reached, continue for some extra number of steps
        if therm_flag
            t_stop += 1
            if t_stop ≥ t_after_therm
                if print_flag
                    println("Stopped at t=$(t)")
                    flush(stdout)
                end
                break
            end
        end
        t += 1
    end
    # push!(evals_all, evals)
    # push!(nb2_all, nb2)
end

serialize("bs/data/bs_L_div_n_10/tt_only/bs_L$(L)_n$(n).dat", 
            Dict(
                #  "evals_all" => evals_all,
                #  "nb2_all" => nb2_all,
                 "L" => L,
                 "n" => n,
                 "interaction_range" => interaction_range,
                 "t_meas" => t_meas,
                 "experiments" => experiments,
                 "tt" => tt,
                 "n_waves" => n_waves,
                #  "roll_window" => roll_window,
                 "therm_fraction" => therm_fraction,
                 "t_after_therm" => t_after_therm))

