include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using LinearAlgebra

function nice_print(s::Vector{Int8})
    println(join(replace(s, 0 => '0', 1 => 'L', 2 => 'A')))
end

function nice_print(states::Vector{Vector{Int8}})
    for s in states
        nice_print(s)
    end
end

function A_local_charge(s::Vector{Int8})
    return Float64.(s .== 2)
end

function A_total_charge(s::Vector{Int8})
    return sum(s .== 2)
end

function A_position(s::Vector{Int8})
    return Float64(findfirst(x -> x == 2, s))
end



# 0 = 0
# 1 = L
# 2 = A
dof_dim = 3

# Construct Hamiltonian on 3 sites
L = 3
pbc = false
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
for site in 1:L-2+2*Int(pbc)
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 1 && s[2] == 2 && s[3] ≠ 2, s -> Int8[1,0,2]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[1,0,2], s -> Int8[1,2,0]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[1,0,2], s -> Int8[1,2,1]))

    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 2 && s[2] ≠ 2 && s[3] ≠ 2 && s ≠ Int8[2,1,1], s -> Int8[2,1,1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[2,1,1], s -> Int8[2,0,0]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[2,1,1], s -> Int8[2,0,1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[2,1,1], s -> Int8[2,1,0]))

    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[1,1,0], s -> Int8[0,0,1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[0,0,1], s -> Int8[1,1,0]))
end
println(length(H_terms)÷2, " terms in the Hamiltonian")
H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);

# Construct the update dictionary
states_all, hams_all = explore_full_space(H, 3)
d = Dict{Vector{Int8}, Vector{Vector{Int8}}}()
for states in states_all
    for s in states
        d[s] = copy(states)
    end
end



# System Parameters
L_list = [10:21;]
L = L_list[parse(Int, ARGS[1])]
# L = 12
pbc = true

function update_function!(state::Vector{Int8})
    classical_brickwork_update!(state, s->true, s -> rand(d[s]); interaction_range=3, pbc=pbc)
end

# Simulation Parameters
experiments = 1000
T_list = [3,5,10,30,50,100,300,500,1000,3000,5000,10000]
T = T_list[parse(Int, ARGS[1])]  # measurement time window
# T = 2
t_max = 2500
measurement_function = A_position
println("L = $(L) \n pbc = $(pbc) \n experiments = $(experiments) \n T = $(T) \n t_max = $(t_max) \n measurement_function = $(measurement_function) ")
flush(stdout)

s_init = Int8[fill(0, L÷2); 2; fill(0, L - (L÷2 + 1))]
evals_all = Vector{typeof(measurement_function(s_init))}[]
init_states = Vector{Int8}[]
for exper in 1:experiments
    println("Experiment: $exper")
    flush(stdout)

    # Initial state preparation
    state = copy(s_init)
    push!(init_states, state)
    println("initial state: $state")
    flush(stdout)

    evals = eltype(eltype(evals_all))[]
    for t in 1:t_max
        if mod(t, 500) == 0
            println("t = $t")
            flush(stdout)
        end

        push!(evals, measure_time_average_and_update!(state; num_samples=T, measurement_function=measurement_function, update_function=update_function!))
        # println("$(t): $(state)")
    end
    push!(evals_all, evals)
end


serialize("data/slow_therm/fuel_AL/L$(L)_T$(T).dat", 
            Dict(
                "L" => L,
                "pbc" => pbc,
                "experiments" => experiments,
                "T" => T,
                "t_max" => t_max,
                "measurement_function" => measurement_function,
                "evals_all" => evals_all,
                "init_states" => init_states
                ))
