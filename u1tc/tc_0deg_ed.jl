include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using LinearAlgebra
using KrylovKit


function coord(x::Int, y::Int, basis_site::Int, basis_sites::Int, Lx::Int, Ly::Int)
    return basis_sites*((mymod(y, Ly)-1)Lx + mymod(x,Lx) - 1) + basis_site
end

function nice_print(s::Vector{Int8})
    # Print top edge with numbers
    println("+" * join(["--$(s[coord(x, 1, 1, 2, Lx, Ly)])--+" for x in 1:Lx]))
    # Print middle rows with numbers
    for y in Ly:-1:1
        # Print row with numbers
        println(join(["$(s[coord(x, y, 2, 2, Lx, Ly)])     " for x in 1:Lx]) * "$(s[coord(1, y, 2, 2, Lx, Ly)])") 
        # Print horizontal edge with numbers
        println("+" * join(["--$(s[coord(x, y, 1, 2, Lx, Ly)])--+" for x in 1:Lx]))
    end
end

function nice_print(states::Vector{Vector{Int8}})
    for (i, s) in enumerate(states)
        println("state $(i):")
        nice_print(s)
        println()
    end
end

function star(x::Int, y::Int)
    return [coord(x, y, 1, 2, Lx, Ly), coord(x, y, 2, 2, Lx, Ly), coord(x-1, y, 1, 2, Lx, Ly), coord(x, y-1, 2, 2, Lx, Ly)]
end

function plaquette(x::Int, y::Int)
    return [coord(x, y, 1, 2, Lx, Ly), coord(x, y, 2, 2, Lx, Ly), coord(x+1, y, 2, 2, Lx, Ly), coord(x, y+1, 1, 2, Lx, Ly)]
end

function row(i::Int)
    return [coord(x, i, 1, 2, Lx, Ly) for x in 1:Lx]
end

function column(i::Int)
    return [coord(i, y, 2, 2, Lx, Ly) for y in 1:Ly]
end

function dual_row(i::Int)
    return [coord(x, i, 2, 2, Lx, Ly) for x in 1:Lx]
end

function dual_column(i::Int)
    return [coord(i, y, 1, 2, Lx, Ly) for y in 1:Ly]
end

function Mz(s::Vector{Int8})
    return -sum(s) + length(s)-sum(s)
end

function Bp(s::Vector{Int8})
    return Int8[(-1)^(reduce(⊻, s[plaquette(x,y)])) for x in 1:Lx, y in 1:Ly]
end

function Wx(s::Vector{Int8})
    if any(Bp(s) .== -1)
        println("Warning: Not all Bp = +1")
    end 
    return (-1)^(reduce(⊻, s[row(1)]))
end

function Wy(s::Vector{Int8})
    if any(Bp(s) .== -1)
        println("Warning: Not all Bp = +1")
    end 
    return (-1)^(reduce(⊻, s[column(1)]))
end

function s_init(Lx::Int, Ly::Int; Wx::Int, Wy::Int)
    s = zeros(Int8, 2Lx*Ly)
    if Wx == -1
        s[dual_column(1)] .= 1
    end
    if Wy == -1
        s[dual_row(1)] .= 1
    end
    return s
end


# 0 = spin-up
# 1 = spin-down
dof_dim = 2
Lx = 2
Ly = 2
λs = -1.0
λp = -1.0
wx = +1
wy = +1


# Construct Hamiltonian
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
for x in 1:Lx, y in 1:Ly
    push!(H_terms, (λs, star(x,y), s -> true, s -> 1 .- s))
end
H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);


# Explore Krylov sector
s = s_init(Lx, Ly; Wx=wx, Wy=wy)
println("Initial state:")
nice_print(s)
states, ham = explore_connected_states(s, H, construct_ham=true, verbose=true);
println("$(length(states)) states found")
flush(stdout)


# Add diagonal H_terms (we know that they are all equal, since Bp's are conserved)
coef = 0
for x in 1:Lx, y in 1:Ly
    coef += λp * (-1)^(reduce(⊻, s[plaquette(x,y)]))
end
if coef ≠ 0
    for i in 1:length(states)
        ham[i, i] = Float64(coef)
    end
end
println("Diagonal terms added.")
flush(stdout)
rows, cols, mels = findnz(ham)


eigenenergies, eigenstates, info = KrylovKit.eigsolve(ham, 6, :SR)
# ground_energy = minimum(eigenenergies)
# ground_states = eigenstates[eigenenergies .≈ ground_energy]

serialize("u1tc/data/tc/Lx$(Lx)_Ly$(Ly)/Wx$(Wx(s))_Wy$(Wy(s)).dat", 
            Dict(
                 "Lx" => Lx,
                 "Ly" => Ly,
                 "λs" => λs,
                 "λp" => λp,
                 "Wx" => Wx(s),
                 "Wy" => Wy(s),
                 "states" => states,
                 "eigenenergies" => eigenenergies,
                 "eigenstates" => eigenstates,
                 "info" => info,
                 "ham" => ham
                 ))



# i = 0
# for state in product(fill(Int8[0:1;], 2Lx*Ly)...)
#     s = collect(state)
#     if all(Bp(s) .== +1) && Mz(s) == 0 && Wx(s) == +1 && Wy(s) == +1
#         nice_print(s)
#         println()
#         i += 1
#     end
# end
# println(i)