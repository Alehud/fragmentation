

"""
    move!(s, H)

Given a state `s` and the Hamiltonian `H`, performs a random nearest neighbor move on the connectivity graph. The whole graph need not be known - the function tries to apply random moves from `H.H_terms`
until it finds one that is applicable. Changes `s` in-place.

# Arguments
- `s::Vector{<:Integer}`: state
- `H::Hamiltonian`: Hamiltonian
"""
function move!(s::Vector{<:Integer}, H::Hamiltonian)
    while true
        _, idx, flippable, flip = rand(H.H_terms)
        if flippable(s[idx])
            s[idx] = flip(s[idx])
            break
        end
    end
end


"""
get_return_times(s_init, H; num_return_times, max_time)

Given an initial state `s_init` and the Hamiltonian `H`, the function performs random moves on the connectivity graph and records times at which the state returns to the initial one.

# Arguments
- `s_init::Vector{<:Integer}`: initial state
- `H::Hamiltonian`: Hamiltonian
- `num_returns::Integer`: desired total number of returns (i.e., `length(return_times)`)
- `max_time::Real=Inf`: time limit after which the function stops (no time limit by default)

# Returns
- `return_times`: vector with return times (each return time is counted from the previous return)
"""
function get_return_times(s_init::Vector{<:Integer}, H::Hamiltonian; num_returns::Integer, max_time::Real=Inf)
    t = 0
    tc = 0
    return_times = Int[]
    s = copy(s_init)
    while t ≤ max_time
        move!(s, H)
        t += 1
        tc += 1
        # println(s)
        if s == s_init
            append!(return_times, tc)
            tc = 0
        end
        if length(return_times) ≥ num_returns
            return return_times
        end
    end
    return return_times
end