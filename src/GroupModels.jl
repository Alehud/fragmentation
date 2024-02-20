export construct_dict, calc_group_elem_representation, calc_group_elem_permutation, construct_d_len_iteratively, word_inverse, isreducible

"""
    construct_dict(L, generator_set, calc_group_elem; type_of_g, omit_identity_letters, no_d, no_identity_states)

Constructs two dictionaries: `d_len` and `d`. The keys of both dictionaries denote group elements. 
The values of `d_len` are the numbers of words (=states) realizing this group element. The values of `d` are vectors with the words themselves.
In addition, the function returns a vector of words realizing the identity group element.

# Arguments
- `L::Integer`: system size
- `generator_set::Vector`: vector with generators (given as arbitrary objects, such as e.g. matrices or functions)
- `calc_group_elem::Function`: function `calc_group_elem(state, generator_set)` that calculates a group element given a state (i.e., sequence of generators) 
    and the generator set
- `omit_identity_letters::Bool = false`: if true, only the words without the identity letter are considered
- `no_d::Bool = false`: if true, the dictionary `d` is not constructed (might be desired for memory considerations)
- `no_identity_states::Bool = false`: if true, the vector `identity_states` is not constructed

# Returns
- `d_len`: dictionary with the number of states
- `d`: dictionary with the states themselves
- `identity_states`: vector of states realizing the identity group element
"""
function construct_dict(L::Integer, generator_set::Vector, calc_group_elem::Function; 
                        omit_identity_letters::Bool = false, no_d::Bool = false, no_identity_states::Bool = false)
    dof_dim = length(generator_set)
    e = calc_group_elem(tuple(zeros(Int8, L)...), generator_set)
    type_of_g = typeof(e)

    d_len = Dict{type_of_g, Integer}()
    d = Dict{type_of_g, Vector{Vector{Int8}}}()
    identity_states = Vector{Int8}[]

    for state in product(fill(Int8[omit_identity_letters:(dof_dim-1);], L)...)
        g = calc_group_elem(state, generator_set)

        d_len[g] = get(d_len, g, 0) + 1
        if !no_d
            d[g] = push!(get(d, g, Vector{Int8}[]), collect(state))
        end
        if !no_identity_states
            if g == e
                push!(identity_states, collect(state))
            end
        end
    end
    return d_len, d, identity_states
end


"""
    calc_group_elem_representation(state, generator_set)

Calculate the group element realized by the `state`, if the generators are represented as matrices

# Arguments
- `state::Union{Tuple{Vararg{Int8}}, Vector{Int8}}`: state (word)
- `generator_set::Vector{<:Matrix}`: vector with generators (given as matrices)

# Returns
- `g`: group element (given as a matrix)
"""
function calc_group_elem_representation(state::Union{Tuple{Vararg{Int8}}, Vector{Int8}}, generator_set::Vector{<:Matrix})
    g = I
    for s in state
        g *= generator_set[s+1]
    end
    return g
end


"""
    calc_group_elem_permutation(state, generator_set, bit_lines)

Calculate the group element realized by the `state`, if the generators are represented as functions acting on an integer from 0 to `2^bit_lines - 1` 
(the functions perform a classical computation on the binary representation of the integer). The group element is then given as a permutation
(a vector `Int8[0:2^bit_lines-1;]` permuted in some way).

# Arguments
- `state::Union{Tuple{Vararg{Int8}}, Vector{Int8}}`: state (word)
- `generator_set::Vector{<:Matrix}`: vector with generators (given as matrices)
- `bit_lines::Integer`: number of bits on which the classical computation is performed

# Returns
- `g`: group element (given as a permutation - a vector `Int8[0:2^bit_lines-1;]` permuted in some way)
"""
function calc_group_elem_permutation(state::Union{Tuple{Vararg{Int8}}, Vector{Int8}}, generator_set::Vector{Function}, bit_lines::Integer)
    bits = Int8[0:2^bit_lines-1;]
    for s in state
        bits = map(generator_set[s+1], bits)
    end
    return bits
end


function construct_d_len_iteratively(L::Tuple{Integer, Integer}, d_len_n::Dict, group_product::Function, filename::String, n::Integer=1)
    group_elems_left = keys(d_len_n)
    g_type = typeof(first(group_elems_left))
    for N in L[1]:L[2]
        println(N)
        flush(stdout)
        d_len_right = deserialize(filename * "_L$(N-n).dat")
        d_len = Dict{g_type, BigInt}()
        for g1 in group_elems_left
            for g2 in keys(d_len_right)
                g = group_product(g1, g2)
                d_len[g] = get(d_len, g, 0) + get(d_len_n, g1, 0) * get(d_len_right, g2, 0)
            end
        end
        serialize(filename * "_L$(N).dat", d_len)
    end
end

function construct_d_len_iteratively(L::Integer, d_len_n::Dict, group_product::Function, all_group_elems::Vector)
    group_elems_1 = keys(d_len_n)
    type_of_g = typeof(first(group_elems_1))

    d_len = Dict{type_of_g, Vector{BigInt}}()
    for g in all_group_elems
        if g in group_elems_1
            d_len[g] = BigInt[d_len_n[g]]
        else
            d_len[g] = [BigInt(0)]
        end
    end

    L_current = length(d_len[first(keys(d_len))])
    for g in all_group_elems
        append!(d_len[g], zeros(BigInt, L - L_current))
    end

    for N in (L_current+1):L
        println(N)
        flush(stdout)
        for g1 in group_elems_1
            for g2 in all_group_elems
                g = group_product(g1, g2)
                d_len[g][N] += get(d_len_n, g1, 0) * d_len[g2][N-1]
            end
        end
    end
    return d_len
end


function construct_d_len_iteratively(L::Integer, d_len_n::Dict, group_product::Function, all_group_elems::Vector, d_len::Dict)
    group_elems_1 = keys(d_len_n)
    type_of_g = typeof(first(group_elems_1))
    L_current = length(d_len[first(keys(d_len))])
    for g in all_group_elems
        append!(d_len[g], zeros(BigInt, L - L_current))
    end
    for N in (L_current+1):L
        println(N)
        flush(stdout)
        for g1 in group_elems_1
            for g2 in all_group_elems
                g = group_product(g1, g2)
                d_len[g][N] += get(d_len_n, g1, 0) * d_len[g2][N-1]
            end
        end
    end
    return d_len
end


"""
    word_inverse(word, inverses)

Compute the inverse of a word, given a dictionary of inverses for all generators (letters) of the presentation.

# Arguments
- `word::Vector{Int8}`: state (word)
- `inverses::Dict{Int8, Int8}`: dictionary that defines an inverse for every letter

# Returns
- inverse of the `word`
"""
function word_inverse(word::Vector{Int8}, inverses::Dict{Int8, Int8})
    return map(x -> get(inverses, x, x), word[end:-1:1])
end


"""
    isreducible(word::Vector{Int8}, relation::Vector{Int8})

Given a generating relation (any word that multiplies to the identity), figure out whether the word can be reduced to a shorter word through a single
application of the relation.

# Arguments
- `word::Vector{Int8}`: state (word)
- `relation::Vector{Int8}`: a word that multiplies to the identity (a generating relation)
- `inverses::Dict{Int8, Int8}`: dictionary that defines an inverse for every letter

# Returns
- `reducible::Bool`: true if the word is reducible
"""
function isreducible(word::Vector{Int8}, relation::Vector{Int8}, inverses::Dict{Int8, Int8})
    len = length(relation)÷2 + 1
    relation_inv = word_inverse(relation, inverses)
    for i in 0:length(relation)-1
        rel = [relation[(i+j) % length(relation) + 1] for j in 0:len-1]
        if issubarray(rel, word)
            return true
        end
        rel = [relation_inv[(i+j) % length(relation_inv) + 1] for j in 0:len-1]
        if issubarray(rel, word)
            return true
        end
    end
    return false
end

"""
    isreducible(word::Vector{Int8}, relation::Vector{Int8})

Given a list of generating relations (words that multiply to the identity), figure out whether the word can be reduced to a shorter word through a single
application of any of the relations.

# Arguments
- `word::Vector{Int8}`: state (word)
- `relations::Vector{Vector{Int8}}`: a vector of words that multiply to the identity (generating relations)
- `inverses::Dict{Int8, Int8}`: dictionary that defines an inverse for every letter

# Returns
- `reducible::Bool`: true if the word is reducible
"""
function isreducible(word::Vector{Int8}, relations::Vector{Vector{Int8}}, inverses::Dict{Int8, Int8})
    # Check if the word is reducible through one application of any of the relations
    for relation in relations
        if isreducible(word, relation, inverses)
            return true
        end
    end
    # If not, then consider application of relations that don't change the length of the word, construct the Hamiltonian consisting of such flips
    # and construct the connectivity graph. Then check if any states in this connectivity graph are reducible.
    H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
    for relation in relations[iseven.(length.(relations))]
        len = length(relation)÷2
        relation_inv = word_inverse(relation, inverses)
        for i in 0:length(relation)-1
            rel = [relation[(i+j) % length(relation) + 1] for j in 0:len-1]
            rel_inv = [relation_inv[(i+j) % length(relation_inv) + 1] for j in 0:len-1]
            for site in 1:length(word) - len + 1
                push!(H_terms, (1.0, [site:site+(len-1);], x -> x == rel, x -> word_inverse(rel, inverses)))
                push!(H_terms, (1.0, [site:site+(len-1);], x -> x == rel_inv, x -> word_inverse(rel_inv, inverses)))
            end
        end 
    end
    H = Hamiltonian(length(keys(inverses)), H_terms, check_hermitian=true);
    states, _ = explore_connected_states(word, H; construct_ham=false, verbose=false)
    for state in states
        for relation in relations
            if isreducible(state, relation, inverses)
                return true
            end
        end
    end
    return false
end
