export explore_connected_states, explore_full_space


"""
    explore_connected_states(s_init, H; construct_ham, check_nonzero, verbose, print_step)

Given an initial state `s_init` and the Hamiltonian `H`, the function iteratively explores the Hilbert space that is connected to the initial state,
and returns the Hamiltonian of this connected subspace, as well as its basis states.

# Arguments
- `s_init::Vector{<:Integer}`: initial state
- `H::Hamiltonian`: Hamiltonian
- `construct_ham::Bool=true`: if true, construct the Hamiltonian (as a sparse matrix)
- `check_nonzero::Bool=false`: if true, check for zero matrix elements and delete them after the Hamiltonian is constructed
- `verbose::Bool=true`: prints the progress if true
- `print_step::Integer=100`: print the progress after every `print_step` states are checked

# Returns
- `states`: vector of basis states of the connected subspace of the Hilbert space. In this basis ham is written.
- `ham`: Hamiltonian as a sparse csr matrix (i.e., for each non-zero matrix element, it stores its row number, column nimber and the value of the element)
"""
function explore_connected_states(s_init::Vector{<:Integer}, H::Hamiltonian; construct_ham::Bool=true, check_nonzero::Bool=false, verbose::Bool=true, print_step::Integer=100)
    # Create a collection of states. At first, only s_init is in the collection.
    states = [s_init]
    count = 0
    rows = Int[]
    cols = Int[]
    mels = Complex[]

    # ## DELETE LATER
    # break_flag = false
    # ##

    # We iterate over each state in the collection and search for states connected to it. Append these states to the collection.
    while count < length(states)
        if verbose
            if (count % print_step == 0) || (count == length(states) - 1)
                println("state: $count, found: $(length(states))")
                flush(stdout)
            end
        end

        state = states[count+1]
        for (coef, idx, flippable, flip) in H.H_terms
            if flippable(state[idx])
                state_new = copy(state)
                state_new[idx] = flip(state_new[idx])
                # Check if we already have such a state in our collection
                if !(state_new in states)
                    # If not, then add it to the collection
                    push!(states, state_new)
                    if construct_ham
                        # Add the corresponding matrix element
                        push!(rows, count)
                        push!(cols, length(states) - 1)
                        push!(mels, coef)
                    end
                elseif construct_ham
                    # If yes, then determine the index of this state
                    ind = findfirst(x -> x == state_new, states) - 1
                    # Add the corresponding matrix element
                    push!(rows, count)
                    push!(cols, ind)
                    push!(mels, coef)
                end

                # ## IMPLEMENT TARGET STATE
                # if state_new == fill(Int8(0), length(s_init))
                #     println("Target state found.")
                #     flush(stdout)
                #     break_flag = true
                #     break
                # end
                # ##
            end
        end
        count += 1

        # ## IMPLEMENT TARGET STATE
        # if break_flag
        #     break
        # end
        # ##

    end

    if construct_ham
        # Construct the Hamiltonian matrix
        if length(mels) > 0
            # This sparse matrix construction will sum up coefficients with similar row and col
            # I.e., if row[i] == row[j] and col[i] == col[j], then the matrix element in the sparse matrix will be mel[i] + mel[j]
            ham = sparse(rows.+1, cols.+1, mels)
            if check_nonzero
                # Remove zero elements
                dropzeros!(ham)
            end
        else
            ham = sparse(Int[], Int[], Complex[])
        end
        return states, ham
    else
        return states, nothing
    end
end


"""
    explore_full_space(H, N_sites; construct_ham)

Explore the full Hilbert space of the model on `N_sites` with Hamiltonian `H`.

# Arguments
- `H::Hamiltonian`: Hamiltonian
- `N_sites::Integer`: total number of sites
- `construct_ham::Bool=true`: if true, construct the Hamiltonian (as a sparse matrix)

# Returns
- `states_all`: vectors of vectors of basis states of the connected subspace of the Hilbert space. In this basis ham is written.
- `hams`: vector of Hamiltonians as sparse csr matrices (i.e., for each non-zero matrix element, it stores its row number, column nimber and the value of the element)
"""
function explore_full_space(H::Hamiltonian, N_sites::Integer; construct_ham::Bool=true)
    states_all = Vector{Vector{Int8}}[]
    hams = SparseMatrixCSC{Complex, Int64}[]
    for state_init_tup in product(fill([0:(H.dof_dim-1);], N_sites)...)
        state_init = collect(state_init_tup)
        if !any([state_init in states for states in states_all])
            states, ham = explore_connected_states(state_init, H, construct_ham=construct_ham, verbose=false)
            push!(states_all, states)
            push!(hams, ham)
        end
    end
    return states_all, hams
end


"""
    explore_full_space(H, states_to_explore; construct_ham)

Explore the Hilbert space of the model connected to `states_to_explore` with Hamiltonian `H`.

# Arguments
- `H::Hamiltonian`: Hamiltonian
- `states_to_explore::Vector{Vector{<:Integer}}`: total number of sites
- `construct_ham::Bool=true`: if true, construct the Hamiltonian (as a sparse matrix)

# Returns
- `states_all`: vectors of vectors of basis states of the connected subspace of the Hilbert space. In this basis ham is written.
- `hams`: vector of Hamiltonians as sparse csr matrices (i.e., for each non-zero matrix element, it stores its row number, column nimber and the value of the element)
"""
function explore_full_space(H::Hamiltonian, states_to_explore::Vector{<:Vector{<:Integer}}; construct_ham::Bool=true)
    states_all = Vector{Vector{Int8}}[]
    hams = []
    for state_init in states_to_explore
        if !any([state_init in states for states in states_all])
            states, ham = explore_connected_states(state_init, H, construct_ham=construct_ham, verbose=false)
            push!(states_all, states)
            push!(hams, ham)
        end
    end
    return states_all, hams
end


