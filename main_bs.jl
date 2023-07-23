include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Graphs
using Karnak
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
L = 6
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

# for H_terms in H.H_terms
#     println(H_terms)
# end

s_init = fill(Int8(0), L)
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
legend = ["0" => "_", "1" => "a", "2" => "b", "3" => "A", "4" => "B"]
state_labels = vertex_labels(states, legend)
remove_vacuum!(state_labels, g)
println("$(length(state_labels)) states left after removing vacuum.")


y = [length(a_star(g, 1, v)) for v in vertices(g)]
window_height = 4000
window_width = 4000
pts = vcat([between.(O + (-window_width/2, window_height/2 - d*window_height/maximum(y)), O + (window_width/2, window_height/2 - d*window_height/maximum(y)), range(0, 1, length=count(i -> i == d, y)+2))[2:end-1] 
            for d in 0:maximum(y)]...)[invperm(sortperm(y))]
# Stress(initialpos = [(pt.x, pt.y) for pt in pts], iterations = 100)


@pdf begin
    background("grey10")
    sethue("orange")
    drawgraph(g, layout = spring,
    vertexshapes = (vtx) -> box(O, 10(L+1), 20, :fill),
    vertexlabels = state_labels,
    vertexlabelfontsizes = 20,
    vertexshapesizes = 12,
    vertexlabeltextcolors = colorant"black",
    vertexlabeloffsetdistances = 0,
    edgelines = (k, s, d, f, t) ->
        arrow(f, t, arrow_func(s, d, y),
            linewidth = 3,
            arrowheadlength = 0),
    edgestrokeweights = 4,
    edgecurvature = 30)
end window_height+100 window_width+100 "pics/bs_test"


@pdf begin
    background("grey10")
    sethue("orange")
    drawgraph(g, layout = pts,
    vertexshapesizes = 4,
    edgelines = (k, s, d, f, t) ->
        arrow(f, t, arrow_func(s, d, y),
            linewidth = 1,
            arrowheadlength = 0),
    edgestrokeweights = 2)
end window_height+100 window_width+100 "pics/bs_nolabel_test"