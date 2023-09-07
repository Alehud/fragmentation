include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Plots
using Colors
using LaTeXStrings
using Statistics
using Polynomials


function thermalize!(state::Vector{Int8}; n_layers::Integer, interaction_range::Integer, generator_set::Vector{Matrix{Union{Int64, Rational{Int64}}}}, d::Dict{Matrix{Union{Rational{Int64}, Int64}}, Vector{Vector{Int8}}})
    L = length(state)
    if interaction_range ≤ L
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
        end
    else
        raise(warning("Interaction range larger than the system size."))
    end
end

function measure_b!(state::Vector{Int8}; n_layers::Integer, interaction_range::Integer, generator_set::Vector{Matrix{Union{Int64, Rational{Int64}}}}, d::Dict{Matrix{Union{Rational{Int64}, Int64}}, Vector{Vector{Int8}}})
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
generator_set = Matrix{Union{Int64, Rational{Int64}}}[e, a, b, ā, b̄]

L = 120
n = L ÷ 10
interaction_range = 3

d = deserialize("data/fragile/d_i$(interaction_range).dat")

t_therm = 0
t_meas = 10^5

experiments = 1
tt = Int64[]
exp_vals_all = Vector{Vector{Float64}}[]
exp_vals_sq_av_all = Vector{Float64}[]

for e in 1:experiments
    println("Experiment: $e")
    flush(stdout)

    state = Int8[fill(2, n); 1; fill(4, n); 3; fill(2, n); 3; fill(4, n); 1]
    while length(state) < L
        insert!(state, rand(1:length(state)), 0)
    end
    # state = Int8.(fill(0, L))
    println("initial state: $state")
    flush(stdout)

    exp_vals = Vector{Float64}[]
    t = 1
    while true
        println("t = $t")
        flush(stdout)
        push!(exp_vals, measure_b!(state, n_layers=t_meas, interaction_range=interaction_range, generator_set=generator_set, d=d))
        if t_therm > 0
            thermalize!(state, n_layers=t_therm, interaction_range=interaction_range, generator_set=generator_set, d=d)
        end
        if t > 100
            break
        end
        # if sum(exp_vals[end].^2) ≤ sum(exp_vals[1].^2)/10
        #     break
        # end
        t += 1
    end
    # exp_vals_sqr = (x -> x.^2).(exp_vals)
    # exp_vals_sqr_av = sum.(exp_vals_sqr) / L

    # t = findfirst(x -> x ≤ exp_vals_sqr_av[1]/10, exp_vals_sqr_av)
    # if isnothing(t)
    #     t = n_cycles
    # end
    push!(tt, t*(t_meas + t_therm))
    push!(exp_vals_all, exp_vals)
    # push!(exp_vals_sq_av_all, exp_vals_sqr_av)
end

serialize("data/fragile/bs_L$(L)_n$(n)_$(int2str(t_therm))_$(int2str(t_meas)).dat", 
            Dict("exp_vals_all" => exp_vals_all,
                 "L" => L,
                 "n" => n,
                 "interaction_range" => interaction_range,
                 "t_therm" => t_therm,
                 "t_meas" => t_meas,
                 "experiments" => experiments,
                 "tt" => tt))

