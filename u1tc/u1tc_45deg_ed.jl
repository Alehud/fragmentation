include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using LinearAlgebra
using KrylovKit

function coord(x::Int, y::Int, Lx::Int, Ly::Int, basis_site::Int=1, basis_sites::Int=1)
    return basis_sites*((mymod(y, Ly)-1)Lx + mymod(x,Lx) - 1) + basis_site
end

function nice_print(s::Vector{Int8})
    for y in Ly:-1:1
        println(join(["   $(isodd(x+y) ? "x" : " ")  " for x in 1:Lx]))
        println(join(["$(s[coord(x, y, Lx, Ly)])     " for x in 1:Lx])) 
    end
end

function nice_print(states::Vector{Vector{Int8}})
    for (i, s) in enumerate(states)
        println("state $(i):")
        nice_print(s)
        println()
    end
end


function plaquette(x::Int, y::Int)
    return [coord(x, y, Lx, Ly), coord(x+1, y, Lx, Ly), coord(x, y+1, Lx, Ly), coord(x+1, y+1, Lx, Ly)]
end

function row(i::Int)
    return [coord(x, i, Lx, Ly) for x in 1:Lx]
end

function column(i::Int)
    return [coord(i, y, Lx, Ly) for y in 1:Ly]
end

function Mz(s::Vector{Int8})
    return -sum(s) + length(s)-sum(s)
end

function Bp(s::Vector{Int8})
    return Int8[(-1)^(reduce(⊻, s[plaquette(x+iseven(y),y)])) for y in Ly:-1:1, x in 1:2:Lx]
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
        raise(error("Odd Lx, Ly are incompatible."))
    end
    if iseven(Lx÷2) && iseven(Ly÷2)
        if Wx == Wy == +1
            return repeat(Int8[fill(1, Lx); fill(0, Lx)], Ly ÷ 2)
        elseif Wx == -1 && Wy == +1
            return repeat(Int8[fill(1, Lx-1); 0; fill(0, Lx-1); 1; 0; fill(1, Lx-1); 1; fill(0, Lx-1);], Ly ÷ 4)
        elseif Wx == +1 && Wy == -1
            return Int8[repeat([1,0,0,1], Lx÷4); repeat([1,0], (Lx÷2)*(Ly-2)); repeat([0,1,1,0], Lx÷4);]
        elseif Wx == Wy == -1
            return Int8[fill(1, Lx-1); 0; repeat([fill(0, Lx-2); 1; 0; 0; fill(1, Lx-3); 0; 0; 0; 1; fill(0, Lx-2); 0; 0; fill(1, Lx-3); 0], (Ly-4)÷4); fill(0, Lx-2); 1; 0; 0; fill(1, Lx-3); 0; 0; 0; fill(1, Lx-1);]
        else
            raise(error("Wx and Wy must be ±1."))
        end
    elseif isodd(Lx÷2) && isodd(Ly÷2)
        if Wx == Wy == +1
            s = fill(Int8(0), Lx*Ly)
            for i in 1:(Ly÷2)-1
                s[row(i)] .= 1
            end
            s[coord(2, Ly÷2, Lx, Ly)] = 1
            s[coord(3, Ly÷2, Lx, Ly)] = 1
            for y in Ly÷2+1:Ly-1
                s[coord(2+y-(Ly÷2+1), y, Lx, Ly)] = 1
            end
            for y in Ly÷2+1:Ly-2
                s[coord(4+y-(Ly÷2+1), y, Lx, Ly)] = 1
            end
            s[coord(Ly÷2+1, Ly-1, Lx, Ly)] = 1
            return s
        elseif Wx == -1 && Wy == +1
            return repeat(Int8[1, 0], (Lx*Ly)÷2)
        elseif Wx == +1 && Wy == -1
            return repeat(Int8[fill(1, Lx); fill(0, Lx)], Ly ÷ 2)
        elseif Wx == Wy == -1
            return Int8[repeat([fill(1, Lx-1); 0; fill(0, Lx-1); 1; 0; fill(1, Lx-1); 1; fill(0, Lx-1);], (Ly-2)÷4); fill(1, Lx-1); 0; fill(0, Lx-1); 1;]
        else
            raise(error("Wx and Wy must be ±1."))
        end
    else
        raise(error("Lx, Ly must be either both even, or both odd."))
    end
end


# ID = parse(Int, ARGS[1])


# 0 = spin-up
# 1 = spin-down
dof_dim = 2
Lx = 8
Ly = 8
λs = -1.0
λp = 0.0
wx_list = [+1, +1, -1]
wy_list = [+1, -1, -1]
wx = -1
wy = -1


# Construct Hamiltonian
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]

for x in 1:Lx, y in 1:Ly
    if isodd(x+y)
        push!(H_terms, (λs, plaquette(x,y), s -> sum(s) == 2, s -> 1 .- s))
    end
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
serialize("u1tc/data/u1tc/45deg/Lx$(Lx)_Ly$(Ly)/Mz$(Mz(s_init))_Wx$(Wx(s_init))_Wy$(Wy(s_init)).dat", 
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

# for wx in [-1, 1], wy in [-1, 1]
#     data = deserialize("u1tc/data/u1tc/45deg/Lx6_Ly6/Mz0_Wx$(wx)_Wy$(wy).dat")
#     eigenenergies = data["eigenenergies"]
#     E0 = minimum(eigenenergies)
#     println("Wx: $(wx), Wy: $(wy), E0 = $(E0)")
# end