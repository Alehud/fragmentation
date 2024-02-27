export apply_single_ham_term!, classical_brickwork_update!, measure_time_average_and_update!, get_return_times

"""
apply_single_ham_term!(state, H)

Given a state `state` and the Hamiltonian `H`, applies a single random Hamiltonian term. I.e., performs a random nearest neighbor move on the connectivity graph. 
The whole graph need not be known - the function tries to apply random moves from `H.H_terms` until it finds one that is applicable (or until `max_attempts` tries). 
Changes `state` in-place.

# Arguments
- `state::Vector{Int8}`: state
- `H::Hamiltonian`: Hamiltonian
- `max_attempts::Integer=1`: maximal number of attempts (the function will try at least length(H.H_terms) times)
"""
function apply_single_ham_term!(state::Vector{Int8}, H::Hamiltonian; max_attempts::Integer=1)
    for _ in 1:max(max_attempts, length(H.H_terms))
        _, idx, flippable, flip = rand(H.H_terms)
        if flippable(state[idx])
            state[idx] = flip(state[idx])
            break
        end
    end
end


"""
classical_brickwork_update!(state, flippable, flip; interaction_range, pbc)

Given a state `state`, performs a classical circuit update with gates arranged in a brickwork pattern. Each gate is of `interaction_range` size. 
It checks if `flippable(state) == true`, and if yes, does `flip(state)`. The number of layers in the circuit is equal to `interaction_range`.
`state` is changed in-place.

# Arguments
- `state::Vector{Int8}`: state
- `flippable::Function`: function that checks if a local pattern is flippable
- `flip::Function`: function that flips a local pattern
- `interaction_range::Integer`: size of local updates
- `pbc::Bool`: pbc/obc
"""
function classical_brickwork_update!(state::Vector{Int8}, flippable::Function, flip::Function; interaction_range::Integer, pbc::Bool)
    L = length(state)
    if interaction_range ≤ L
        for layer in 1:interaction_range
            start_site = rand(1:L)
            idx = mymod([start_site:start_site+interaction_range-1;], L)
            for _ in 1:(L - (!pbc)*(layer-1))÷interaction_range
                if flippable(state[idx])
                    state[idx] = flip(state[idx])
                end
                if pbc
                    idx = mymod(idx .+ interaction_range, L)
                else
                    idx = idx .+ interaction_range
                end
            end
        end
    else
        raise(error("Interaction range larger than the system size."))
    end
end



"""
measure_time_average_and_update!(state; num_samples, measurement_function, update_function!, num_updates_between_measurements=1)

Performs a time-averaged measurement over `num_samples` samples using `measurement_function`. The updates are performed using `update_function`, with `num_updates_between_measurements` updates between two consecutive collected samples.

# Arguments
- `state::Vector{Int8}`: state
- `num_samples::Integer`: number of measurement samples
- `measurement_function::Function`: function used to perform a measurement
- `update_function!::Function`: function used to perform an update (must change a state in-place)
- `num_updates_between_measurements::Integer=1`: number of updates between two measurements
"""
function measure_time_average_and_update!(state::Vector{Int8}; num_samples::Integer, measurement_function::Function, update_function::Function, num_updates_between_measurements::Integer=1)
    observable_sum = measurement_function(state)
    for _ in 1:num_updates_between_measurements
        update_function(state)
    end
    for _ in 1:(num_samples-1)
        for _ in 1:num_updates_between_measurements
            update_function(state)
        end
        observable_sum += measurement_function(state)
    end
    return observable_sum / num_samples
end





"""
get_return_times(s_init, H; num_return_times, max_time)

Given an initial state `s_init` and the Hamiltonian `H`, the function performs random moves on the connectivity graph and records times at which the state returns to the initial one.

# Arguments
- `s_init::Vector{<:Integer}`: initial state
- `H::Hamiltonian`: Hamiltonian
- `move_function!::Function`: function use to perform a single move
- `num_returns::Integer`: desired total number of returns (i.e., `length(return_times)`)
- `max_time::Real=Inf`: time limit after which the function stops (no time limit by default)

# Returns
- `return_times`: vector with return times (each return time is counted from the previous return)
"""
function get_return_times(s_init::Vector{<:Integer}, H::Hamiltonian, move_function!::Function; num_returns::Integer, max_time::Real=Inf)
    t = 0
    tc = 0
    return_times = Int[]
    s = copy(s_init)
    while t ≤ max_time
        move_function!(s, H)
        t += 1
        tc += 1
        # println(s)
        if s == s_init
            append!(return_times, tc)
            tc = 0
            r = length(return_times)
            if r % 1000 == 0
                println("$r returns") 
                flush(stdout)
            end
            if r ≥ num_returns
                return return_times
            end
        end
    end
    return return_times
end



