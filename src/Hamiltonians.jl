export Hamiltonian


"""
Hamiltonian
# Fields
- `dof_dim::Integer`: each d.o.f. is represented by an integer from 0 to `dof_dim`-1
- `H_terms::Vector{Tuple}`: terms of the Hamiltonian, vector of tuples of the form `(coef, idx, flippable, flip)`, where `coef` is a coefficient (number),
                            `idx` is a vector of indices of sites (integers, counting starts from 1), `flippable` and `flip` are functions that take a vector of
                            length `length(idx)` as an argument; `flippable` returns a boolean variable that determines whether the local configuration is flippable; 
                            `flip` is a function that returns a vector with a new configuration.

"""
struct Hamiltonian
    dof_dim::Integer
    H_terms::Vector{Tuple}
    
    function Hamiltonian(dof_dim::Integer, H_terms::Vector{Tuple{Real, Vector{<:Integer}, Function, Function}}; check_hermitian::Bool = false, reduce_hamiltonian::Bool = false)
        for (coef, idx, flippable, flip) in H_terms
            if !(typeof(coef) <: Number)
                raise(error("`coef` should be a real or complex number."))
            end
            if !(typeof(idx) <: Vector{<:Integer})
                raise(error("`idx` should be a vector of integers."))
            end
            if !(flippable isa Function)
                raise(error("`flippable` should be a function."))
            end
            if !(flip isa Function)
                raise(error("`flip` should be a function."))
            end
            if !(typeof(flippable(fill(0, length(idx)))) <: Bool)
                raise(error("`flippable` should return a boolean."))
            end
            if !(typeof(flip(fill(0, length(idx)))) <: Vector{<:Integer})
                raise(error("`flip` should return a vector of length `length(idx)` with integers."))
            end
            for state in [collect(x) for x in product(fill([0:(dof_dim-1);], length(idx))...)]
                if flippable(state)
                    state_new = flip(state)
                    if any(state_new .≥ dof_dim) || any(state_new .< 0)
                        raise(error("`flip` function should not return values outside of [0, `dof_dim`-1] range when applied to flippable configurations."))
                    end
                end
            end
        end
        if reduce_hamiltonian
            deleteat!(H_terms, [term[1] == 0.0 for term in H_terms])
            if length(H_terms) > 1
                i = 1
                j = 2
                while true
                    coef, idx, flippable, flip = H_terms[i]
                    coef1, idx1, flippable1, flip1 = H_terms[j]
                    if idx == idx1 && flippable == flippable1 && flip == flip1
                        if coef + coef1 ≠ 0.0
                            H_terms[i] = (coef + coef1, idx, flippable, flip)
                            deleteat!(H_terms, j)
                        else
                            deleteat!(H_terms, j)
                            deleteat!(H_terms, i)
                            j = i + 1
                        end
                    else
                        j += 1
                    end
                    if j > length(H_terms)
                        i += 1
                        j = i + 1
                    end
                    if i ≥ length(H_terms)
                        break
                    end
                end
            end
        end
        if check_hermitian
            for (coef, idx, flippable, flip) in H_terms
                for state in [collect(x) for x in product(fill([0:(dof_dim-1);], length(idx))...)]
                    if flippable(state)
                        flag = true
                        state_new = flip(state)
                        for (coef1, idx1, flippable1, flip1) in H_terms
                            if idx == idx1 && flippable1(state_new) && flip1(state_new) == state && coef1 == conj(coef)
                                flag = false
                                break
                            end
                        end
                        if flag
                            raise(error("Hamiltonian is not Hermitian."))
                        end
                    end
                end
            end
        end
        new(dof_dim, H_terms)
    end
end