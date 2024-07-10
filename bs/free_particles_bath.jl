include("../src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Statistics
# using StatsBase


# 0 = 0
# 1 = b

L = 100
n = L
interaction_range = 3
T = 20
experiments = 10000
pbc = false

println(" L = $(L) \n n = $n \n T = $(T) \n experiments = $(experiments) \n")
flush(stdout)

# Define update function
function update_function!(state::Vector{Int8})
    state[1] = Int8(rand(0:1))  # bath on the first site flips the state randomly
    classical_brickwork_update!(state, s -> sum(s) == 1, s -> Int8(1) .- s; interaction_range=2, pbc=pbc)
end
# Define measurement function
function nb(state::Vector{Int8})
    return state .== 1
end


tt = Int64[]    # thermalization times
evals_all = Vector{Vector{Float64}}[]
nb2_all = Vector{Float64}[]

# roll_window = 100
therm_fraction = 0.1    # if nb2 reaches this fraction of the initial value, then thermalization occurs
t_after_therm = 100     # if thermalization occured, continue for this time and then stop
t_max = 2500    # stop after this time no matter what
print_rate = experiments÷1000 + Int(experiments < 1000)
for exper in 1:experiments
    print_flag = mod(exper, print_rate) == 0
    if print_flag
        println("----------------------------------------------------------------------------")
        println("Experiment: $exper")
        flush(stdout)
    end

    # Initial state preparation
    # all-b word (b^n) sprinkled with identities 
    state = Int8[fill(1, n);]
    while length(state) < L
        insert!(state, rand(1:length(state)), 0)
    end

    # init_state = copy(state)
    if print_flag
        println("initial state: $state")
        flush(stdout)
    end

    # evals = eltype(eltype(evals_all))[]
    nb2 = eltype(eltype(nb2_all))[]
    
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

        eval = measure_time_average_and_update!(state; num_samples=T, measurement_function=nb, update_function=update_function!)
        # push!(evals, eval)
        push!(nb2, sum(eval.^2)/L)

        # Stop after maximal time
        if t > t_max
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
                push!(tt, t*T)
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
        #         push!(tt, t*T)
        #         therm_flag = true
        #         println("Thermalized at t=$(t) (with roll_window=$(roll_window))")
        #         flush(stdout)
        #     end
        # end
        # If thermalization occured, continue for some extra number of time steps
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
    push!(nb2_all, nb2)
end

serialize("bs/data/free_particles_bath/free_particles_bath_L$(L)_n$(n)_T$(int2str(T)).dat", 
            Dict(
                #  "evals_all" => evals_all,
                 "nb2_all" => nb2_all,
                 "L" => L,
                 "n" => n,
                 "interaction_range" => interaction_range,
                 "T" => T,
                 "experiments" => experiments,
                 "tt" => tt,
                #  "roll_window" => roll_window,
                 "therm_fraction" => therm_fraction,
                 "t_after_therm" => t_after_therm))
        