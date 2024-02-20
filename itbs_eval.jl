include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Statistics


function thermalize!(state::Vector{Int8}; n_layers::Integer, interaction_range::Integer, d::Dict{Vector{Int8}, Vector{Vector{Int8}}})
    L = length(state)
    if interaction_range ≤ L
        for layer in 1:n_layers
            for sublayer in 1:interaction_range
                idx = [sublayer:sublayer+interaction_range-1;]
                while idx[end] ≤ L
                    state[idx] = rand(d[state[idx]])
                    idx = idx .+ interaction_range
                end
            end
        end
    else
        raise(warning("Interaction range larger than the system size."))
    end
end

function measure_c!(state::Vector{Int8}; n_layers::Integer, interaction_range::Integer, d::Dict{Vector{Int8}, Vector{Vector{Int8}}})
    L = length(state)
    if interaction_range ≤ L
        state_sum = fill(Int64(0), L)
        for layer in 1:n_layers
            for sublayer in 1:interaction_range
                idx = [sublayer:sublayer+interaction_range-1;]
                while idx[end] ≤ L
                    state[idx] = rand(d[state[idx]])
                    idx = idx .+ interaction_range
                end
            end
            state_sum = state_sum .+ (state .== 5) .- (state .== 6)
        end
    else
        raise(warning("Interaction range larger than the system size."))
    end
    return state_sum / n_layers
end


# 0 = e
# 1 = a
# 2 = b
# 3 = a^(-1)
# 4 = b^(-1)
# 5 = c
# 6 = c^(-1)


L_list = [66:69;]
# n_list = [12,10,9,8,6]
# t_meas_list = [10^5, 10^4, 3*10^3, 10^3, 10^2]


L = L_list[parse(Int, ARGS[1])]
n = 4
n_waves = 1
interaction_range = 3

d = deserialize("data/itbs/d_L$(interaction_range).dat")

t_therm = 0
t_meas = 100

experiments = 1000
tt = Int64[]
exp_vals_all = Vector{Vector{Float64}}[]
exp_vals_sq_av_all = Vector{Float64}[]

println(" L = $(L) \n n = $n \n n_waves = $(n_waves) \n t_meas = $(t_meas) \n experiments = $(experiments) \n")

for exper in 1:experiments
    println("Experiment: $exper")
    flush(stdout)

    # Construct the initial state
    w = Int8[fill(5, n); 4; fill(6, n);]
    w_inv = Int8[fill(5, n); 2; fill(6, n);]
    wp = Int8[w; 1; w_inv;]
    wp_inv = Int8[w; 3; w_inv;]
    state = repeat(Int8[wp_inv; 1; wp; 3;], n_waves)
    # state = fill(Int8(0), L)
    while length(state) < L
        insert!(state, rand(1:length(state)), 0)
    end
    println("initial state: $state")
    flush(stdout)

    print("t = 0 ")
    println(count(x -> x==5, state), " ", count(x -> x==6, state))
    flush(stdout)

    exp_vals = Vector{Float64}[]
    t = 1
    t_stop = 0
    therm_flag = false
    while true
        if mod(t, 1000) == 0
            print("t = $t ")
            println(count(x -> x==5, state), " ", count(x -> x==6, state))
            flush(stdout)
        end
        push!(exp_vals, measure_c!(state, n_layers=t_meas, interaction_range=interaction_range, d=d))
        if t_therm > 0
            thermalize!(state, n_layers=t_therm, interaction_range=interaction_range, d=d)
        end
        if t > 100000
            push!(tt, -1)
            break
        end
        if therm_flag == false
            if count(x -> x==5, state) == 0
                push!(tt, t*(t_meas + t_therm))
                therm_flag = true
            end
        end
        if therm_flag
            t_stop += 1
        end
        if t_stop ≥ 100
            break
        end
        t += 1
    end
    # exp_vals_sqr = (x -> x.^2).(exp_vals)
    # exp_vals_sqr_av = sum.(exp_vals_sqr) / L

    # t = findfirst(x -> x ≤ exp_vals_sqr_av[1]/10, exp_vals_sqr_av)
    # if isnothing(t)
    #     t = n_cycles
    # end
    # push!(exp_vals_all, exp_vals)
end

serialize("data/itbs/itbs_n$(n)/only_tt/itbs_L$(L)_n$(n)_$(int2str(t_meas)).dat", 
            Dict(
                #  "exp_vals_all" => exp_vals_all,
                 "L" => L,
                 "n" => n,
                 "interaction_range" => interaction_range,
                 "t_therm" => t_therm,
                 "t_meas" => t_meas,
                 "experiments" => experiments,
                 "tt" => tt,
                 "n_waves" => n_waves))