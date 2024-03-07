include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization


# 0 = e
# 1 = a
# 2 = b
# 3 = a^(-1)
# 4 = b^(-1)
# 5 = c
# 6 = c^(-1)
dof_dim = 7
e = Int64[1 0; 0 1]
a = Int64[1 1; 0 1]
ā = Int64[1 -1; 0 1]
b = Int64[2 0; 0 1]
b̄ = [1//2 0; 0 1]
generator_set_bs = Matrix{Rational{Int64}}[e, a, b, ā, b̄]

d_bs_L1 = deserialize("data/bs/d_L1.dat")
d_bs_L2 = deserialize("data/bs/d_L2.dat")
d_bs_L3 = deserialize("data/bs/d_L3.dat")

# Construct Hamiltonian
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
for (group_element, states) in d_bs_L3
    len = length(states)
    for i in 1:len, j in 1:len
        if i ≠ j
            push!(H_terms, (1.0, [1,2,3], x -> x == states[i], x -> states[j]))
        end
    end
    push!(H_terms, (1.0, [1,2,3], x -> x == Int8[2,2,5], x -> Int8[5,2,0]))
    push!(H_terms, (1.0, [1,2,3], x -> x == Int8[5,2,0], x -> Int8[2,2,5]))
    push!(H_terms, (1.0, [1,2,3], x -> x == Int8[4,4,5], x -> Int8[5,4,0]))
    push!(H_terms, (1.0, [1,2,3], x -> x == Int8[5,4,0], x -> Int8[4,4,5]))
    push!(H_terms, (1.0, [1,2,3], x -> x == Int8[6,2,2], x -> Int8[2,6,0]))
    push!(H_terms, (1.0, [1,2,3], x -> x == Int8[2,6,0], x -> Int8[6,2,2]))
    push!(H_terms, (1.0, [1,2,3], x -> x == Int8[6,4,4], x -> Int8[4,6,0]))
    push!(H_terms, (1.0, [1,2,3], x -> x == Int8[4,6,0], x -> Int8[6,4,4]))
end
for site in 1:2
    for (group_element, states) in d_bs_L2
        len = length(states)
        for i in 1:len, j in 1:len
            if i ≠ j
                push!(H_terms, (1.0, [site, site+1], x -> x == states[i], x -> states[j]))
            end
        end
    end
    push!(H_terms, (1.0, [site, site+1], s -> s==Int8[5,0] || s==Int8[6,0] || s==Int8[0,5] || s==Int8[0,6], s -> [s[2],s[1]]))
    push!(H_terms, (1.0, [site, site+1], s -> s==Int8[5,6] || s==Int8[6,5], s -> Int8[0,0]))
end
H = Hamiltonian(dof_dim, H_terms, check_hermitian=false);

# Construct dictionary
d = Dict{Vector{Int8}, Vector{Vector{Int8}}}()
for state in product(fill([0:6;], 3)...)
    states, _ = explore_connected_states(Int8.(collect(state)), H; construct_ham=false, check_nonzero=false)
    d[Int8.(collect(state))] = states
end


# Do some checks
for (key, value) in d
    if any(key .== 5) || any(key .== 6)
        println(key, " -> ", value)
    end
    if all(key .< 5)
        g = reduce(*, generator_set_bs[key .+ 1])
        for s in value
            if !(s in d_bs_L3[g])
                println("aaaaa ", key)
            end
        end
        for s in d_bs_L3[g]
            if !(s in value)
                println("AAAAA ", key)
            end
        end
    end
end


serialize("data/itbs/d_L3.dat", d)

