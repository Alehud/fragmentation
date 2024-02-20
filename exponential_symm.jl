include("src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using LinearAlgebra
# using Graphs
# using Luxor
# using Karnak
# using Colors
# using Plots

function charge(s::Vector{Int8})
    x = [0:length(s)-1;]
    mz = replace(s, 0 => -1, 1 => 0, 2 => 1)
    return sum(@. 2^x * mz)
end

function nice_print(s::Vector{Int8})
    println(join(replace(s, 0 => '-', 1 => '_', 2 => '+')))
end

function nice_print(states::Vector{Vector{Int8}})
    for s in states
        print("Q=", charge(s), "  ")
        nice_print(s)
    end
end

function semipositive_state(nz::Vector{<:Integer}, trail_zeros::Integer=0)
    s = Int8[]
    if trail_zeros > 0
        append!(s, Int8[fill(1, trail_zeros);])
    end
    for n in nz
        append!(s, Int8[2; fill(1, n);])
    end
    return s
end

function alternating_state(nz::Vector{<:Integer}, trail_zeros::Integer=0)
    s = Int8[]
    if trail_zeros > 0
        append!(s, Int8[fill(1, trail_zeros);])
    end
    i = 1
    for n in nz
        append!(s, Int8[1+i; fill(1, n);])
        i *= -1
    end
    return s
end

function num_states_T_semipositive(nz::Vector{<:Integer}, pbc::Bool)
    if pbc
        num_states = [1 0 0; 0 1 0; 0 0 1;]
    else
        num_states = [1, nz[1], 1]
    end
    for n in nz[2-Int(pbc):end-1+Int(pbc)]
        T = [1 1 0
             n n n
             1 1 1]
        num_states = T * num_states
    end
    if pbc
        num_states = tr(num_states)
    else
        n = nz[end]
        T = [1 1 0
             n n n
             0 0 0]
        num_states = T * num_states
        num_states = ([1 1 1] * num_states)[1]
    end
    return num_states
end

function num_states_T_alternating(nz::Vector{<:Integer}, pbc::Bool)
    if isodd(length(nz)) && pbc
        println("NEED TO MODIFY TRANSFER MATRIX FOR ODD NZ WITH PBC")
        return nothing
    end
    if pbc
        num_states = [1 0 0; 0 1 0; 0 0 1;]
    else
        num_states = [1, nz[1], 1]
    end
    for n in nz[2-Int(pbc):end-1+Int(pbc)]
        T = [1 1 1
             n n 0
             1 1 0]
        num_states = T * num_states
    end
    if pbc
        num_states = tr(num_states)
    else
        n = nz[end]
        T = [1 1 1
             n n 0
             0 0 0]
        num_states = T * num_states
        num_states = ([1 1 1] * num_states)[1]
    end
    return num_states
end



# (S_i^-)^2 S_{i+1}^+  +  (S_i^+)^2 S_{i+1}^-
# 0 - spin-down
# 1 - spin-0
# 2 - spin-up

dof_dim = 3
pbc = true

nz = [2, 3, 4, 3]
s_init = alternating_state(nz)
println("s_init: $(s_init)")
println("Transfer matrix approach: $(num_states_T_alternating(nz, pbc))")

println(prod(nz) + nz[1]*nz[2] + nz[2]*nz[3] + nz[3]*nz[4] + nz[4]*nz[1] + 1 + 1)


# Construct Hamiltonian
L = length(s_init)
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
for site in 1:L-1+Int(pbc)
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1], x -> (x[1] == 0) && (x[2] > 0), x -> [x[1] + 2, x[2] - 1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1], x -> (x[1] == 2) && (x[2] < 2), x -> [x[1] - 2, x[2] + 1]))
end
println(length(H_terms)÷2, " terms in the Hamiltonian")
H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);

# Explore Krylov sector
states, ham = explore_connected_states(s_init, H, construct_ham=true, verbose=false);
println("$(length(states)) states found")
# rows, cols, mels = findnz(ham)

nice_print(states)

if isreal(ham)
    ham = real.(ham)
else
    throw(error("Hamiltonian is not real."))
end
println("Hamiltonian: ")
display(ham)
ham = Matrix(ham)
energies = eigvals(ham)
vecs = eigvecs(ham)
println("Eigenenergies: ", round.(energies, digits=4))
println("Eigenvectors: ")
for v in eachcol(vecs)
    println(round.(v, digits=4))
end



for s_init in product(fill(Int8[0,1,2], L)...)
    states, ham = explore_connected_states(collect(s_init), H, construct_ham=true, verbose=false);
    flag = false
    states_trim = (s -> s[s .≠ 1]).(states)
    n = count((s -> all(s .≠ circshift(s, 1))).(states_trim))
    if n ≠ 1
        println(n)
        nice_print(states)
    end
    # if !flag
    #     println(s_init)
    # end
end