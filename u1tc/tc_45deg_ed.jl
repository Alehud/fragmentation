include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using LinearAlgebra
using KrylovKit

ID = parse(Int, ARGS[1])

function job()

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

    function s_init(Lx::Int, Ly::Int; Wx::Int, Wy::Int)
        if isodd(Lx) || isodd(Ly)
            raise(error("Odd Lx, Ly are incompatible."))
        end
        s = zeros(Int8, Lx*Ly)
        if Wx == -1
            s[column(1)] = 1 .- s[column(1)] 
        end
        if Wy == -1
            s[row(1)] = 1 .- s[row(1)]
        end
        return s
    end


    # 0 = spin-up
    # 1 = spin-down
    dof_dim = 2
    Lx = 6
    Ly = 6
    λs = -1.0
    λp = -1.0
    wx_list = [+1, +1, -1, -1]
    wy_list = [+1, -1, +1, -1]
    wx = wx_list[ID]
    wy = wy_list[ID]


    # Construct Hamiltonian
    H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]

    for x in 1:Lx, y in 1:Ly
        if isodd(x+y)
            push!(H_terms, (λs, plaquette(x,y), s -> true, s -> 1 .- s))
        end
    end
    H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);


    # Explore Krylov sector
    s = s_init(Lx, Ly; Wx=wx, Wy=wy)
    println("Initial state:")
    nice_print(s)
    println("Mz: $(Mz(s)) \nWx: $(Wx(s)) \nWy: $(Wy(s)) \nBp: $(Bp(s))")
    states, ham = explore_connected_states(s, H, construct_ham=true, verbose=true);
    println("$(length(states)) states found")
    flush(stdout)

    # Add diagonal H_terms (we know that they are all equal, since Bp's are conserved)
    coef = 0
    for x in 1:Lx, y in 1:Ly
        if isodd(x+y)
            coef += λp * (-1)^(reduce(⊻, s[plaquette(x,y)]))
        end
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
    serialize("u1tc/data/tc/45deg/Lx$(Lx)_Ly$(Ly)/Wx$(Wx(s))_Wy$(Wy(s)).dat", 
                Dict(
                    "Lx" => Lx,
                    "Ly" => Ly,
                    "λs" => λs,
                    "λp" => λp,
                    "Mz" => Mz(s),
                    "Wx" => Wx(s),
                    "Wy" => Wy(s),
                    "states" => states,
                    "eigenenergies" => eigenenergies,
                    "eigenstates" => eigenstates,
                    "info" => info,
                    "ham" => ham
                    ))
end

job()

# for wx in [-1, 1], wy in [-1, 1]
#     data = deserialize("u1tc/data/tc/45deg/Lx6_Ly6/Wx$(wx)_Wy$(wy).dat")
#     eigenenergies = data["eigenenergies"]
#     E0 = minimum(eigenenergies)
#     println("Wx: $(wx), Wy: $(wy), E0 = $(E0)")
# end