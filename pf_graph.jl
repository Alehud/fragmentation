include("src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using LinearAlgebra
using Graphs
using Luxor
using Karnak
using Colors
using Statistics
using Plots


function get_krylov_label(state::Vector{Int8})
    label = copy(state)  # Create a copy of the input state to avoid modifying it
    i = 1
    while i < length(label)
        if label[i] == label[i+1]
            splice!(label, i:i+1)  # Remove adjacent equal elements
            i = max(1, i - 1)  # Move one step back to check for more adjacent equal elements
        else
            i += 1
        end
    end
    return label
end

function convert_to_points(points::Vector; radius::Real)
    pts_max = sqrt(maximum((p -> p[1]^2 + p[2]^2).(points)))
    return [Point(p[1], p[2])*radius/(2*pts_max) for p in points]
end

# function label_to_pos(label::Vector{Int8}, L::Integer)
#     if isempty(label)
#         return Point(0.0, 0.0)
#     end
#     layer = (length(label) + isodd(L)) ÷ 2
#     nr = count(c -> c == 0, label)
#     ng = count(c -> c == 1, label)
#     nb = count(c -> c == 2, label)

# end


# 0 = red
# 1 = green
# 2 = blue
dof_dim = 3

# Construct Hamiltonian
L = 5
H_terms = Tuple[]
for site in 1:(L-1)
    for i in 1:dof_dim-1
        push!(H_terms, (1.0, [site, site+1], x -> x[1] == x[2], x -> mod.(x .+ i, dof_dim)))
    end
end
for i in 1:dof_dim-1
    push!(H_terms, (1.0, [1], x -> true, x -> mod.(x .+ i, dof_dim)))
end
H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);


s_init = fill(Int8(0), L)
println("s_init: $s_init")
states, ham = explore_connected_states(s_init, H, construct_ham=true);
println("$(length(states)) states found", )
rows, cols, mels = findnz(ham)
krylov_labels = get_krylov_label.(states)
different_krylov_labels = unique(krylov_labels)
krylov_adjacency_matrix = [label1 == label2 for label1 in krylov_labels, label2 in krylov_labels]
krylov_adjacency_matrix[diagind(krylov_adjacency_matrix)] .= false
krylov_adjacency_matrix = SparseMatrixCSC{Int64, Int64}(krylov_adjacency_matrix)


edges = Edge.(collect(zip(rows, cols)))
g = SimpleGraph(edges)
legend = ["0" => "r", "1" => "g", "2" => "b"]
state_labels = vertex_labels(states, legend)

# adj_orange = krylov_adjacency_matrix .* adjacency_matrix(g)
# adj_red = (1 .- krylov_adjacency_matrix) .* adjacency_matrix(g)
# g_orange = SimpleGraph(adj_orange)
# g_red = SimpleGraph(adj_red)
# distance_matrix_orange = [length(a_star(g_orange, v1, v2)) for v1 in 1:nv(g_orange), v2 in 1:nv(g_orange)]
# distance_matrix_red = [length(a_star(g_red, v1, v2)) for v1 in 1:nv(g_red), v2 in 1:nv(g_red)]
# distance_matrix_between = [length(a_star(g, v1, v2))*isempty(a_star(g_orange, v1, v2)) for v1 in 1:nv(g), v2 in 1:nv(g)]
# weights_orange = ifelse.(distance_matrix_orange .== 0, 0, distance_matrix_orange.^(-2))
# weights_red = ifelse.(distance_matrix_red .== 0, 0, distance_matrix_red.^(-2))
# weights_between = ifelse.(distance_matrix_between .== 0, 0, distance_matrix_between.^(-2))


window_height = 4000
window_width = 4000            
pts = convert_to_points(stress(adjacency_matrix(g)), radius=min(window_height, window_width))


for kr_lab in different_krylov_labels
    vertices = findall(lab -> lab == kr_lab, krylov_labels)
    radius = length(vertices)*30
    vertices_pos = copy(pts[vertices])
    average_pos = mean(vertices_pos)
    sg, vmap = induced_subgraph(g, vertices)
    if length(vertices) > 1
        pts_sg = convert_to_points(stress(adjacency_matrix(sg)), radius=radius) .+ average_pos
    else
        pts_sg = [Point(0.0, 0.0)] .+ average_pos
    end
    pts[vertices] = pts_sg
end

@pdf begin
    background("grey10")
    sethue("orange")
    drawgraph(g, layout = pts,
    vertexshapes = (vtx) -> ellipse(O, 10(L+2), 20, :fill),
    vertexshapesizes = 12,
    vertexlabels = state_labels,
    vertexlabelfontsizes = 20,
    vertexlabeltextcolors = colorant"black",
    vertexlabeloffsetdistances = 0,
    vertexfillcolors = (vtx) -> distinguishable_colors(length(different_krylov_labels), colorant"black"; dropseed=true)[findfirst(x -> x == krylov_labels[vtx], different_krylov_labels)],
    edgestrokeweights = 4,
    edgestrokecolors = (edgenumber, to, from, src, dest) -> begin 
                                                                if krylov_labels[to] ≠ krylov_labels[from] 
                                                                    return colorant"red"
                                                                else
                                                                    return colorant"orange"
                                                                end
                                                            end
    )
end window_height+100 window_width+100 "pics/pf/test_L$(L)"



