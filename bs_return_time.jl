include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Graphs
using Karnak
using Plots
using Colors

# 0 = vacuum
# 1 = a
# 2 = b
# 3 = a^(-1)
# 4 = b^(-1)

function annihilatable(s)
    return s[1]*s[2] == 3 || s[1]*s[2] == 8
end

function annihilate(s)
    return [0, 0]
end

function creatable(s)
    return s[1] == 0 && s[2] == 0
end

function create_a1(s)
    return [1, 3]
end

function create_a2(s)
    return [3, 1]
end

function create_b1(s)
    return [2, 4]
end

function create_b2(s)
    return [4, 2]
end

function flippable(s)
    return s[1] == 2 && s[2] == 1 && s[3] == 0
end

function flip(s)
    return [1,1,2]
end

function flippable_back(s)
    return s[1] == 1 && s[2] == 1 && s[3] == 2
end

function flip_back(s)
    return [2,1,0]
end

function hoppable(s)
    return s[1] * s[2] == 0 && s[1] + s[2] ≠ 0
end

function hop(s)
    return [s[2], s[1]]
end


dof_dim = 5
L = 10
H_terms = Tuple[]
for i in 1:(L-2)
    push!(H_terms, (1.0, [i, i+1, i+2], flippable, flip))
    push!(H_terms, (1.0, [i, i+1, i+2], flippable_back, flip_back))
end
for i in 1:(L-1)
    push!(H_terms, (1.0, [i, i+1], annihilatable, annihilate))
    push!(H_terms, (1.0, [i, i+1], creatable, create_a1))
    push!(H_terms, (1.0, [i, i+1], creatable, create_a2))
    push!(H_terms, (1.0, [i, i+1], creatable, create_b1))
    push!(H_terms, (1.0, [i, i+1], creatable, create_b2))
    push!(H_terms, (1.0, [i, i+1], hoppable, hop))
end

H = Hamiltonian(dof_dim, H_terms, check_hermitian=true)


s_init = fill(Int8(0), L)
println("s_init: $s_init")

num_return_times = 10
return_times = get_return_times(s_init, H, move!; num_returns=num_return_times)

serialize("data/fragile/bs_return_times_L$(L)_$(num_return_times ÷ 1000)k.dat", return_times)

# histogram(return_times)
# savefig("pics/bs_return_times_$(num_return_times)")

# return_times = deserialize("data/fragile/bs_return_times_L$(L)_$(num_return_times ÷ 1000)k.dat")