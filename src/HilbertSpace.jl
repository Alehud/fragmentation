export explore_connected_states


"""
    explore_connected_states(s_init, H_terms, N_sites, check_hermitian, return_ham)

Given an initial state `s_init` and the Hamiltonian terms `H_terms`, the function iteratively explores the Hilbert space that is connected to the initial state,
and returns the Hamiltonian of this connected subspace, as well as its basis states.

# Arguments
- `s_init::Integer`: initial state
- `H_terms::Vector{Tuple}`: the list of terms in the Hamiltonian
- `N_sites::Integer`: the number of sites in the system
- `construct_ham::Bool=true`: if true, construct the Hamiltonian as a sparse matrix
- `check_nonzero::Bool=true`: if true, check for zero matrix elements and delete them after the Hamiltonian is constructed
- `check_hermitian::Bool=false`: check if the Hamiltonian is Hermitian, raise a warning if not

# Returns
- `ham`: Hamiltonian as a sparse csr matrix (i.e., for each non-zero matrix element, it stores its row number, column nimber and the value of the element)
- `states`: list of basis states of the connected subspace of the Hilbert space. In this basis ham is written.
- `rows`: list of row numbers for all non-zero matrix elements of ham
- `cols`: list of column numbers for all non-zero matrix elements of ham
- `mels`: list of values of all non-zero matrix elements
"""
function explore_connected_states(s_init::Integer, H_terms::Vector{Tuple}, N_sites::Integer; construct_ham::Bool=true, check_nonzero::Bool=true, check_hermitian::Bool=false)
    # Total number of bits in a state (TODO: add more than 2 d.o.f. at each site)
    N_tot = N_sites

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
        end

        state = states[count+1]
        for term in H_terms
            coef, ops = term
            s_temp = state
            annihilation_flag = false
            for op in reverse(ops)
                site, dag = op
                if get_bit(s_temp, site) ⊻ (dag == '+')  # check if the state is annihilated
                    s_temp = flip_bit(s_temp, site)  # Switch 0->1 or 1->0
                else
                    annihilation_flag = true
                    break
                end
            end
            if !annihilation_flag
                # Check if we already have such a state in our collection
                if !(s_temp in states)
                    # If not, then add it to the collection
                    push!(states, s_temp)
                    if construct_ham
                        # Add the corresponding matrix element
                        push!(rows, count)
                        push!(cols, length(states) - 1)
                        push!(mels, coef)
                    end
                elseif construct_ham
                    # If yes, then determine the index of this state
                    ind = findfirst(x -> x == s_temp, states) - 1
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
                # Access new nonzero matrix elements and their row and column indices
                rows, cols, mels = findnz(ham)
                rows .-= 1
                cols .-= 1
            end
        else
            ham = sparse(Int[], Int[], Complex[])
        end
    
        if check_hermitian
            if nnz(ham - adjoint(ham)) > 0
                throw(Warning("Hamiltonian is not Hermitian!"))
            end
        end
    else
        ham = nothing
        rows = nothing
        cols = nothing
        mels = nothing
    end
    return states, ham, rows, cols, mels
end
