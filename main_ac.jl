include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Graphs
using Karnak
using Colors


function flippable(s)
    return s[1] == 2 && s[2] == 1 && s[3] == 0
end

function flip(s)
    return [1, 2, 2]
end

function flippable_back(s)
    return s[1] == 1 && s[2] == 2 && s[3] == 2
end

function flip_back(s)
    return [2, 1, 0]
end

function hoppable(s)
    return s[1] * s[2] == 0 && s[1] + s[2] ≠ 0
end

function hop(s)
    return [s[2], s[1]]
end


dof_dim = 3
L = 12
H_terms = Tuple[]
for i in 1:(L-2)
    push!(H_terms, (1.0, [i, i+1, i+2], flippable, flip))
    push!(H_terms, (1.0, [i, i+1, i+2], flippable_back, flip_back))
end
for i in 1:(L-1)
    push!(H_terms, (1.0, [i, i+1], hoppable, hop))
end

H = Hamiltonian(dof_dim, H_terms, check_hermitian=true)

# for H_terms in H.H_terms
#     println(H_terms)
# end

s_init = Int8[2,1,1,1,0,0,0,0,0,0,0,0]
println("s_init: $s_init")
states, ham = explore_connected_states(s_init, H, construct_ham=true)
println("$(length(states)) states found", )
# println("states:")
# for state in states
#     println(state)
# end
# println("ham: ", ham)
display(ham)
rows, cols, mels = findnz(ham)
# println("rows: $rows")
# println("cols: $cols")
# println("mels: $mels")

# states_all, hams = explore_full_space(H, 8, construct_ham=true)
# for (i, states) in enumerate(states_all)
#     println(i, " ", states)
# end

# rows = []
# cols = []
# mels = []

# function f()
#     i = 0
#     for (ham, states) in zip(hams, states_all)
#         rows_temp, cols_temp, mels_temp = findnz(ham)
#         append!(rows, rows_temp .+ i)
#         append!(cols, cols_temp .+ i)
#         append!(mels, mels_temp)
#         i += length(states)
#     end
#     print(i)
#     println(rows)
#     println(cols)
# end

# f()


edges = Edge.(collect(zip(rows, cols)))
g = SimpleGraph(edges)
legend = ["0" => "_", "1" => "a", "2" => "c"]
state_labels = vertex_labels(states, legend)
remove_vacuum!(state_labels, g)
println("$(length(state_labels)) states left after removing vacuum.")

y = [length(a_star(g, 1, v)) for v in vertices(g)]
H = 2000
W = 2000
pts = vcat([between.(O + (-W/2, H/2 - d*H/maximum(y)), O + (W/2, H/2 - d*H/maximum(y)), range(0, 1, length=count(i -> i == d, y)+2))[2:end-1] 
            for d in 0:maximum(y)]...)[invperm(sortperm(y))]
# Stress(initialpos = [(pt.x, pt.y) for pt in pts], iterations = 100)

function arrow_func(s, d)
    if abs(y[s]-y[d]) == 1 || abs(s-d) == 1
        return [0, 0]
    else
        return [20, 25]
    end
end

@pdf begin
    background("grey10")
    sethue("orange")
    drawgraph(g, layout = pts,
    vertexshapes = (vtx) -> box(O, 10(L+1), 20, :fill),
    vertexlabels = state_labels,
    vertexlabelfontsizes = 20,
    vertexshapesizes = 12,
    vertexlabeltextcolors = colorant"black",
    vertexlabeloffsetdistances = 0,
    edgelines = (k, s, d, f, t) ->
        arrow(f, t, arrow_func(s, d),
            linewidth = 3,
            arrowheadlength = 0),
    edgestrokeweights = 4,
    edgecurvature = 30)
end H+100 W+100 "pics/ac"


# Mz = (N_sites - sum(s_init_arr))*1 + sum(s_init_arr)*(-1)
# println("Mz: ", Mz)
# serialize("data/Lx=$(Lx)_Ly=$(Ly)_Mz=$(Mz).dat", (length(states), states))
# (N_states, states) = deserialize("Lx=$(Lx)_Ly=$(Ly)_Mz=0.dat")
# println(N_states)
# brint(states, N_sites)



