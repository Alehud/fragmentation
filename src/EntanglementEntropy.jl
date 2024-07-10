export reduced_density_matrix, bipartite_entanglement_entropy, entanglement_entropy_of_eigenstates_in_connected_sector

"""
    reduced_density_matrix(rho, states, bipartition_mask)

Calculates the reduced density matrix of density matrix `rho` in region A encoded as `true` in `bipartition_mask`. `states` contain the product state basis.

# Arguments
- `rho::Union{Matrix{<:Number}, SparseMatrixCSC{<:Number, <:Integer}}`: density matrix
- `states::Vector{Vector{Int8}}`: vector with product states comprising the basis
- `bipartition_mask::Vector{Bool}`: bipartition mask (`true` denotes region A, `false` denotes region B)

# Returns
- `rho_a`: reduced density matrix
"""
function reduced_density_matrix(rho::Union{Matrix{<:Number}, SparseMatrixCSC{<:Number, <:Integer}}, states::Vector{Vector{Int8}}, bipartition_mask::Vector{Bool})
    # D.o.f. of each state in `states` in A region (can contain equal elements). It is a vector with |phi_i^A>.
    states_a = (s -> s[bipartition_mask]).(states)
    # D.o.f. of each state in `states` in B region (can contain equal elements). It is a vector with |phi_i^B>.
    states_b = (s -> s[.!bipartition_mask]).(states)

    # |phi_i^A> might be equal to |phi_j^A> for i not equal to j. We have to take care of this and sum up rows and columns
    # corresponding to equal |phi_i^A>.
    # idx_lists_states_a is a list of lists of indices corresponding to equal states |phi_i^A>.
    # For example, [[1,3], [2], [4,5,6]] would mean that |phi_1^A> = |phi_3^A>, |phi_4^A> = |phi_5^A> = |phi_6^A>.
    idx_lists_states_a = [findall(x -> x == state, states_a) for state in unique(states_a)]

    if typeof(rho) <: Matrix    # if rho is a normal Matrix
        # Boolean matrix showing which elements in `states_b`` are equal to each other, basically it is <phi_i^B|phi_j^B>
        states_b_eq_map = [s1 == s2 for s1 in states_b, s2 in states_b]
        # rho_a is a reduced density matrix, rho_a = \sum_{i,j} c_i c_j^* B_{ij} |phi_i^A><phi_j^A|
        rho_a = rho .* states_b_eq_map
        # Sum rows and columns corresponding to equal |phi_i^A>.
        rho_a = vcat([sum(rho_a[idxs, :], dims=1) for idxs in idx_lists_states_a]...)
        rho_a = hcat([sum(rho_a[:, idxs], dims=2) for idxs in idx_lists_states_a]...)

    else    # if rho is a sparse matrix
        # Decompose sparse matrix into indices
        rows, cols, mels = findnz(rho)

        # Set the coefficient to 0 if <phi_i^B|phi_j^B> = 0
        idx_to_leave = [states_b[r] == states_b[c] for (r,c) in zip(rows, cols)]
        rows = rows[idx_to_leave]
        cols = cols[idx_to_leave]
        mels = mels[idx_to_leave]

        # Sum rows and columns corresponding to equal |phi_i^A>.
        for idx in idx_lists_states_a
            replace!(x -> x in idx ? idx[1] : x, rows);
            replace!(x -> x in idx ? idx[1] : x, cols);
            (i -> i - sum(i .> idx[2:end])).(rows)
            (i -> i - sum(i .> idx[2:end])).(cols)
        end
        # Compose the sparse matrix back
        rho_a = sparse(cols, rows, mels)
    end
    return rho_a
end


"""
    bipartite_entanglement_entropy(rho, states, bipartition_mask, threshold_zero_eigvals=1e-10)

Computes bipartite entanglement entropy for the density matrix `rho`. This function is inefficient for pure states -- instead, use method `bipartite_entanglement_entropy(psi, states, bipartition_mask, threshold_zero_eigvals)` 
    or `bipartite_entanglement_entropy(Ψ, threshold_zero_eigvals)`.

# Arguments
- `rho::Union{Matrix{<:Number}, SparseMatrixCSC{<:Number, <:Integer}}`: density matrix
- `states::Vector{Vector{Int8}}`: vector with product states comprising the basis
- `bipartition_mask::Vector{Bool}`: bipartition mask (`true` denotes region A, `false` denotes region B)
- `threshold_zero_eigvals::Float64=1e-10`: eigenvalues of rho_a with absolute value below this threshold are ignored 

# Returns
- `entropy`: bipartite entanglement entropy
"""
function bipartite_entanglement_entropy(rho::Union{Matrix{<:Number}, SparseMatrixCSC{<:Number, <:Integer}}, states::Vector{Vector{Int8}}, bipartition_mask::Vector{Bool}, threshold_zero_eigvals::Float64=1e-10)
    rho_a = reduced_density_matrix(rho, states, bipartition_mask)
    rho_a_eigvals = eigvals(Matrix(rho_a))

    # Get rid of 0 eigenvalues (eigenvalues of a density matrix are always non-negative)
    rho_a_eigvals = rho_a_eigvals[rho_a_eigvals .> threshold_zero_eigvals]
    entropy = - sum(rho_a_eigvals .* log.(rho_a_eigvals))
    return entropy
end


"""
    bipartite_entanglement_entropy(Ψ, threshold_zero_eigvals=1e-10)

Computes bipartite entanglement entropy for a pure state with density matrix rho = |psi><psi|. State |psi> given as a matrix Ψ, which has basis states from region A as rows and basis states from region B as columns.

# Arguments
- `Ψ::Matrix{<:Number}`: matrix, which has basis states from region A as rows and basis states from region B as columns (the list of basis states need not span the full Hilbert space)
- `threshold_zero_eigvals::Float64=1e-10`: eigenvalues of the reduced density matrix with absolute value below this threshold are ignored 

# Returns
- `entropy`: bipartite entanglement entropy
"""
function bipartite_entanglement_entropy(Ψ::Matrix{<:Number}, threshold_zero_eigvals::Float64=1e-10)
    # Find singular values of Ψ. Their squares are eigenvalues of the reduced density matrix
    # rho_a_eigvals = (svdvals(Ψ)).^2
    rho_a_eigvals = abs2.(LinearAlgebra.LAPACK.gesvd!('N', 'N', Ψ)[2])
    

    # Get rid of 0 eigenvalues (eigenvalues of a density matrix are always non-negative)
    rho_a_eigvals = rho_a_eigvals[rho_a_eigvals .> threshold_zero_eigvals]
    entropy = - sum(rho_a_eigvals .* log.(rho_a_eigvals))
    return entropy
end


"""
    bipartite_entanglement_entropy(psi, states, bipartition_mask, threshold_zero_eigvals=1e-10)

Computes bipartite entanglement entropy for a pure state with density matrix rho = |psi><psi|.

# Arguments
- `psi::Vector{<:Number}`: state
- `states::Vector{Vector{Int8}}`: vector with product states comprising the basis
- `bipartition_mask::Vector{Bool}`: bipartition mask (`true` denotes region A, `false` denotes region B)
- `threshold_zero_eigvals::Float64=1e-10`: eigenvalues of the reduced density matrix with absolute value below this threshold are ignored 

# Returns
- `entropy`: bipartite entanglement entropy
"""
function bipartite_entanglement_entropy(psi::Vector{<:Number}, states::Vector{Vector{Int8}}, bipartition_mask::Vector{Bool}, threshold_zero_eigvals::Float64=1e-10)
    # D.o.f. of each state in `states` in A region (can contain equal elements). It is a vector with |phi_i^A>.
    states_a = (s -> s[bipartition_mask]).(states)
    # D.o.f. of each state in `states` in B region (can contain equal elements). It is a vector with |phi_i^B>.
    states_b = (s -> s[.!bipartition_mask]).(states)

    # Bases for A and B regions
    basis_a = unique(states_a)
    basis_b = unique(states_b)

    # Construct matrix Ψ, where rows/columns denote product states in A/B region.
    Ψ = zeros(typeof(psi[1]), length(basis_a), length(basis_b))
    for (c, sa, sb) in zip(psi, states_a, states_b)
        row = findfirst(s -> s == sa, basis_a)
        col = findfirst(s -> s == sb, basis_b)
        Ψ[row, col] = c
    end

    return bipartite_entanglement_entropy(Ψ, threshold_zero_eigvals)
end


"""
    entanglement_entropy_of_eigenstates_in_connected_sector(ham, states, bipartition_mask)

Computes entanglement entropy of eigenstates in a dynamically connected sector.

# Arguments
- `ham::SparseMatrixCSC{T, Ti}`: Hamiltonian as a sparse matrix
- `states::Vector{Vector{Int8}}`: vector with product states comprising the basis
- `bipartition_mask::Vector{Bool}`: bipartition mask (`true` denotes region A, `false` denotes region B)

# Returns
- `entropies`: bipartite entanglement entropy
- `eigvals`: eigenenergies
- `eigvecs`: eigenstates
"""
function entanglement_entropy_of_eigenstates_in_connected_sector(ham::SparseMatrixCSC{T, Ti}, states::Vector{Vector{Int8}}, bipartition_mask::Vector{Bool}) where {T, Ti<:Integer}
    eigvals, eigvecs = eigen(Matrix(ham))
    entropies = Float64[]

    for psi in eachcol(eigvecs)
        entropy = bipartite_entanglement_entropy(Vector(psi), states, bipartition_mask)
        push!(entropies, entropy)
    end

    return entropies, eigvals, eigvecs
end

