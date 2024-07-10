include("../src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Statistics
# using StatsBase


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

L_list = [80:5:150;]
T_list = [Int(round(16*1.25^i)) for i in 0:length(L_list)-1]
experiments_list = [30, 25, 23, 19, 15, 15, 13, 11, 9, 8, 7, 7, 6, 5, 4] * 10^3

L = L_list[ID]
n = L
interaction_range = 3
T = T_list[ID]
experiments = experiments_list[ID]
pbc = false

println(" L = $(L) \n n = $n \n T = $(T) \n experiments = $(experiments) \n")
flush(stdout)

# Download update dictionary
d = deserialize("bs/data/d/d_L$(interaction_range).dat")

# Define update function
function update_function!(state::Vector{Int8})
    state[1] = Int8(rand(0:4))  # bath on the first site flips the state randomly
    classical_brickwork_update!(state, s->true, s -> rand(d[calc_group_elem_representation(s, generator_set)]); interaction_range=3, pbc=pbc)
end
# Define measurement function
function nb(state::Vector{Int8})
    return (state .== 2) .- (state .== 4)
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
    # w_large_half word (b^n a b^{-n}) sprinkled with identities 
    # state = Int8[fill(2, n); 1; fill(4, n);]
    # while length(state) < L
    #     insert!(state, rand(1:length(state)), 0)
    # end
    # all-b word (b^n) sprinkled with identities 
    state = Int8[fill(2, n);]
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
    # push!(nb2_all, nb2)
end

serialize("bs/data/bs_bath/all_b/tt_only/bs_L$(L)_n$(n)_T$(int2str(T)).dat", 
            Dict(
                #  "evals_all" => evals_all,
                #  "nb2_all" => nb2_all,
                 "L" => L,
                 "n" => n,
                 "interaction_range" => interaction_range,
                 "T" => T,
                 "experiments" => experiments,
                 "tt" => tt,
                #  "roll_window" => roll_window,
                 "therm_fraction" => therm_fraction,
                 "t_after_therm" => t_after_therm))


