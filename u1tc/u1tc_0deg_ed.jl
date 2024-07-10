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

function s_init_Mz0(Lx::Int, Ly::Int; Wx::Int, Wy::Int)
    if isodd(Lx) || isodd(Ly)
        raise(error("Only even Lx, Ly are supported."))
    end
    if Wx == Wy == +1
        return repeat(Int8[0, 1], Lx*Ly)
    elseif Wx == -1 && Wy == +1
        s = fill(Int8(0), 2Lx*Ly)
        for x in 2:Lx
            s[column(x)] .= 1
        end
        for y in 1:Ly
            if isodd(y)
                s[coord(1, y, 1, 2, Lx, Ly)] = 1
            else
                s[coord(Lx, y, 1, 2, Lx, Ly)] = 1
            end
        end
        return s
    elseif Wx == +1 && Wy == -1
        s = fill(Int8(0), 2Lx*Ly)
        for y in 2:Ly
            s[row(y)] .= 1
        end
        for x in 1:Lx
            if isodd(x)
                s[coord(x, Ly, 2, 2, Lx, Ly)] = 1
            else
                s[coord(x, 1, 2, 2, Lx, Ly)] = 1
            end
        end
        return s
    elseif Wx == Wy == -1
        s = fill(Int8(0), 2Lx*Ly)
        for y in 1:Ly-1
            s[dual_row(y)] .= 1
        end
        s[dual_column(1)] .= 1
        return s
    else
        raise(error("Wx and Wy must be ±1."))
    end
end


# ID = parse(Int, ARGS[1])


# 0 = spin-up
# 1 = spin-down
dof_dim = 2
Lx = 4
Ly = 4
λs = -1.0
λp = 0.0
# wx_list = [+1, +1, -1]
# wy_list = [+1, -1, -1]
wx = +1
wy = +1
    

# Construct Hamiltonian
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
for x in 1:Lx, y in 1:Ly
    push!(H_terms, (λs, star(x,y), s -> sum(s) == 2, s -> 1 .- s))
end
H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);


# Explore Krylov sector
s_init = s_init_Mz0(Lx, Ly; Wx=wx, Wy=wy)
println("Initial state:")
nice_print(s_init)
println("Mz: $(Mz(s_init)) \nWx: $(Wx(s_init)) \nWy: $(Wy(s_init)) \nBp: $(Bp(s_init))")
states, ham = explore_connected_states(s_init, H, construct_ham=true, verbose=true);
println("$(length(states)) states found")
flush(stdout)

# Add diagonal H_terms (we know that they are all equal, since Bp's are conserved)
# coef = 0
# for x in 1:Lx, y in 1:Ly
#     coef += λp * (-1)^(reduce(⊻, s_init[plaquette(x,y)]))
# end
# if coef ≠ 0
#     for i in 1:length(states)
#         ham[i, i] = Float64(coef)
#     end
# end
# println("Diagonal terms added.")
# flush(stdout)
# rows, cols, mels = findnz(ham)


eigenenergies, eigenstates, info = KrylovKit.eigsolve(ham, 6, :SR)

serialize("u1tc/data/u1tc/0deg/Lx$(Lx)_Ly$(Ly)/Mz$(Mz(s_init))_Wx$(Wx(s_init))_Wy$(Wy(s_init)).dat", 
            Dict(
                "Lx" => Lx,
                "Ly" => Ly,
                "λs" => λs,
                "λp" => λp,
                "Mz" => Mz(s_init),
                "Wx" => Wx(s_init),
                "Wy" => Wy(s_init),
                "states" => states,
                "eigenenergies" => eigenenergies,
                "eigenstates" => eigenstates,
                "info" => info,
                "ham" => ham
                ))

