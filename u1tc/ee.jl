include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using LinearAlgebra
using KrylovKit
using Optimization, OptimizationOptimJL, OptimizationBBO

# model = "tc"
# angle = 0
# Lx = Ly = 2
# gs_degeneracy = 4
# w_list = [(-1,-1), (-1,1), (1,-1), (1,1)]
# w_list = [(-1,-1), (1,1)]
# w_list = [(-1,1), (1,-1), (1,1)]
# w_list = [(-1,-1), (-1,1), (1,-1)]



function job(model::String, angle::Int, Lx::Int, Ly::Int, gs_degeneracy::Int, w_list::Vector{Tuple{Int, Int}})
    if angle == 0
        bipartition_mask_x = [fill(true, Lx*Ly); fill(false, Lx*Ly)]
        bipartition_mask_y = repeat([fill(true, Lx); fill(false, Lx)], Ly)
    elseif angle == 45
        bipartition_mask_x = [fill(true, Lx*Ly ÷ 2); fill(false, Lx*Ly ÷ 2)]
        bipartition_mask_y = repeat([fill(true, Lx÷2); fill(false, Lx÷2)], Ly)
    end

    if model == "tc"
        mz = ""
    elseif model == "u1tc"
        mz = "Mz0_"
    end

    print("Deserializing states...   ")
    flush(stdout)
    gs_states = [Float64[] for i in 1:gs_degeneracy]
    gs_energies = Float64[]
    basis_states = Vector{Int8}[]
    for (i, (wx, wy)) in enumerate(w_list)
        
        data = deserialize("u1tc/data/$(model)/$(angle)deg/Lx$(Lx)_Ly$(Ly)/$(mz)Wx$(wx)_Wy$(wy).dat")
        states = data["states"]
        eigenenergies = data["eigenenergies"]
        gs_energy = minimum(eigenenergies)
        gs_state = data["eigenstates"][eigenenergies .≈ gs_energy]
        if length(gs_state) ≠ 1
            @warn "There should be only 1 ground state in each sector"
        end
        gs_state = gs_state[1]
        push!(gs_energies, gs_energy)
        append!(gs_states[i], gs_state)
        for j in 1:gs_degeneracy
            if j ≠ i
                append!(gs_states[j], zeros(typeof(gs_states[1][1]), length(states)))
            end
        end
        append!(basis_states, states)
    end
    println("States deserialized.")
    println("GS energies: ", gs_energies)
    flush(stdout)


    if model == "tc"
        # Analytical results
        s1x = 1/sqrt(2) * (gs_states[4] + gs_states[3])
        s2x = 1/sqrt(2) * (gs_states[4] - gs_states[3])
        s3x = 1/sqrt(2) * (gs_states[2] + gs_states[1])
        s4x = 1/sqrt(2) * (gs_states[2] - gs_states[1])

        s1y = 1/sqrt(2) * (gs_states[4] + gs_states[2])
        s2y = 1/sqrt(2) * (gs_states[4] - gs_states[2])
        s3y = 1/sqrt(2) * (gs_states[3] + gs_states[1])
        s4y = 1/sqrt(2) * (gs_states[3] - gs_states[1])

        println("Analytical S-matrix: ")
        S_anal = [round(dot(sx, sy); digits=3) for sx in [s1x,s2x,s3x,s4x], sy in [s1y,s2y,s3y,s4y]] * 2
        println(S_anal)
        flush(stdout)
    end



    function spherical_to_euclidean(ϕ::Vector{Float64})
        return [[cos(ϕ[i])*prod(sin, ϕ[1:i-1], init=1) for i in 1:length(ϕ)]; prod(sin, ϕ, init=1);]
    end

    # Define the objective function to minimize
    function objective_function(ϕ::Vector{<:Number}, params)
        gs_basis_matrix, states_a, states_b, basis_a, basis_b = params
        c = spherical_to_euclidean(ϕ)
        psi = gs_basis_matrix * c   # wavefunction
        # Construct matrix Ψ, where rows/columns denote product states in A/B region.
        rows = [findfirst(isequal(s), basis_a) for s in states_a]
        cols = [findfirst(isequal(s), basis_b) for s in states_b]
        Ψ = Matrix(sparse(rows, cols, psi))
        return bipartite_entanglement_entropy(Ψ)
    end

    # Define the optimization function
    # optf = OptimizationFunction(objective_function, Optimization.AutoZygote())
    optf = OptimizationFunction(objective_function)

    function MES(gs_degeneracy::Integer, gs_states::Vector{Vector{Float64}}, states::Vector{Vector{Int8}}, bipartition_mask::Vector{Bool})
        mes = Vector{Float64}[]
        # Matrix with gs states as rows
        gs_basis_matrix = hcat(gs_states...)
        # D.o.f. of each state in `states` in A region (can contain equal elements). It is a vector with |phi_i^A>.
        states_a = (s -> s[bipartition_mask]).(states)
        # D.o.f. of each state in `states` in B region (can contain equal elements). It is a vector with |phi_i^B>.
        states_b = (s -> s[.!bipartition_mask]).(states)
        # Bases for A and B regions
        basis_a = unique(states_a)
        basis_b = unique(states_b)
        for n in gs_degeneracy:-1:2
            println("n: $n")
            ϕ0 = [rand(n-2)*pi; rand(1)*2pi;]
            params = [gs_basis_matrix, states_a, states_b, basis_a, basis_b]

            # Define the optimization problem
            prob = OptimizationProblem(optf, ϕ0, params, lb=fill(0.0, n-1) , ub=[fill(pi, n-2); 2pi;], sense=MinSense)
            # Solve the optimization problem
            # sol = solve(prob, BFGS(); reltol=1e-3)
            sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); abstol=1e-8)
            println(sol.stats)

            c = spherical_to_euclidean(sol.u)
            s = gs_basis_matrix * c
            push!(mes, s)
            gs_basis_matrix[:, 1] = s
            Q, _ = qr(gs_basis_matrix)
            gs_basis_matrix = Matrix(Q)[:, 2:end]
        end
        push!(mes, gs_basis_matrix[:, 1])

        return mes
    end

    println("Computing MES_x...")
    flush(stdout)
    mes_x = MES(gs_degeneracy, gs_states, basis_states, bipartition_mask_x);
    println("Computing MES_y...")
    flush(stdout)
    mes_y = MES(gs_degeneracy, gs_states, basis_states, bipartition_mask_y);
    println("MES computed.")
    flush(stdout)

    if model == "tc"
        println("Comparing to analytical results: ")
        println([isapprox(mes, s; atol=1e-3) || isapprox(mes, -s; atol=1e-3) for mes in mes_x, s in [s1x, s2x, s3x, s4x]])
        println([isapprox(mes, s; atol=1e-3) || isapprox(mes, -s; atol=1e-3) for mes in mes_y, s in [s1y, s2y, s3y, s4y]])
        flush(stdout)
    end

    println("S-matrix (not normalized by the total quantum dimension): ")
    S = [round(dot(sx, sy); digits=3) for sx in mes_x, sy in mes_y]
    println(S)
    flush(stdout)

    print("Serializing results...   ")
    flush(stdout)
    path = "u1tc/data/$(model)/$(angle)deg/Lx$(Lx)_Ly$(Ly)"
    serialize(joinpath(mkpath(path), "ee.dat"), 
                Dict(
                    "Lx" => Lx,
                    "Ly" => Ly,
                    "angle" => angle,
                    "gs_degeneracy" => gs_degeneracy,
                    "bipartition_mask_x" => bipartition_mask_x,
                    "bipartition_mask_y" => bipartition_mask_y,
                    "mes_x" => mes_x,
                    "mes_y" => mes_y,
                    "S" => S
                    ))
    println("Results serialized.")
    flush(stdout)

end


ID = parse(Int, ARGS[1])

model_list = ["u1tc", "u1tc"]
angle_list = [0, 45]
L_list = [4, 6]
gs_degeneracy_list = [2, 3]
w_list_list = [[(-1,-1), (1,1)], 
               [(-1,-1), (-1,1), (1,-1)]]

job(model_list[ID], angle_list[ID], L_list[ID], L_list[ID], gs_degeneracy_list[ID], w_list_list[ID])


