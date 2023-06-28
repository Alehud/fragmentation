export explore_connected_states, explore_full_space


"""
    explore_connected_states(s_init, H_terms, N_sites, check_hermitian, return_ham)

Given an initial state `s_init` and the Hamiltonian terms `H_terms`, the function iteratively explores the Hilbert space that is connected to the initial state,
and returns the Hamiltonian of this connected subspace, as well as its basis states.

# Arguments
- `s_init::Vector{<:Integer}`: initial state
- `H::Hamiltonian`: Hamiltonian
- `construct_ham::Bool=true`: if true, construct the Hamiltonian as a sparse matrix
- `check_nonzero::Bool=true`: if true, check for zero matrix elements and delete them after the Hamiltonian is constructed

# Returns
- `ham`: Hamiltonian as a sparse csr matrix (i.e., for each non-zero matrix element, it stores its row number, column nimber and the value of the element)
- `states`: list of basis states of the connected subspace of the Hilbert space. In this basis ham is written.
- `rows`: list of row numbers for all non-zero matrix elements of ham
- `cols`: list of column numbers for all non-zero matrix elements of ham
- `mels`: list of values of all non-zero matrix elements
"""
function explore_connected_states(s_init::Vector{<:Integer}, H::Hamiltonian; construct_ham::Bool=true, check_nonzero::Bool=false)
    # Create a collection of states. At first, only s_init is in the collection.
    states = [s_init]
    count = 0
    rows = Int[]
    cols = Int[]
    mels = Complex[]

    # We iterate over each state in the collection and search for states connected to it. Append these states to the collection.
    while count < length(states)
        if (count % 100 == 0) || (count == length(states) - 1)
            println("state: $count, found: $(length(states))")
            flush(stdout)
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
            end
        end
        count += 1
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


function explore_full_space(H::Hamiltonian, N_sites::Integer; construct_ham::Bool=true)
    states_all = Vector{Vector{<:Integer}}[]
    hams = []
    for state_init in [collect(x) for x in product(fill([0:(H.dof_dim-1);], N_sites)...)]
        if !any([state_init in states for states in states_all])
            states, ham = explore_connected_states(state_init, H, construct_ham=construct_ham)
            push!(states_all, states)
            push!(hams, ham)
        end
    end
    return states_all, hams
end


