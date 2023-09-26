include("src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using Graphs
using Luxor
using Karnak
using Colors
using Plots

# 0 = vacuum
# 1 = a
# 2 = b
# 3 = a^(-1)
# 4 = b^(-1)

interaction_range = 3
dof_dim = 5

d = deserialize("data/fragile/d_i$(interaction_range).dat")

# Construct Hamiltonian
L = 10
H_terms = Tuple[]
for site in 1:(L-(interaction_range-1))
    for (group_element, states) in d
        len = length(states)
        for i in 1:len, j in 1:len
            if i ≠ j
                push!(H_terms, (1.0, [site:site+(interaction_range-1);], x -> x == states[i], x -> states[j]))
            end
        end
    end
end
println(length(H_terms)÷2, " terms in the Hamiltonian")
H = Hamiltonian(dof_dim, H_terms, check_hermitian=false);


n = 1
s_init = Int8[fill(2, n); 1; fill(4, n); 3; fill(2, n); 3; fill(4, n); 1; fill(0, L - (4n+4))]
println("s_init: $s_init")
states, ham = explore_connected_states(s_init, H, construct_ham=true);
println("$(length(states)) states found", )
rows, cols, mels = findnz(ham)

edges = Edge.(collect(zip(rows, cols)))
g = SimpleGraph(edges)
legend = ["0" => "_", "1" => "a", "2" => "b", "3" => "A", "4" => "B"]
state_labels = vertex_labels(states, legend)
remove_vacuum!(state_labels, g)
println("$(length(state_labels)) states left after removing vacuum.")

target = ""
v_target = findfirst(x -> x == target, state_labels)
path = a_star(g, 1, v_target)
println("Shortest path: ")
for edge in path
    print(state_labels[edge.src], "  ->  ")
end
println(state_labels[v_target])


# d = deserialize("data/fragile/d_i$(L).dat")
# gg = SimpleGraph{Int64}[]
# legend = ["0" => "_", "1" => "a", "2" => "b", "3" => "ā", "4" => "b̄"]
# for (perm, states_symm_sector) in d
#     # Connectivity within symmetry sector
#     states_symm_sector_list, hams = explore_full_space(H, states_symm_sector; construct_ham=true)
#     g = SimpleGraph()
#     i = 0
#     for (states, ham) in zip(states_symm_sector_list, hams)
#         rows, cols, mels = findnz(ham);
#         rows .+= i
#         cols .+= i
#         edges = Edge.(collect(zip(rows, cols)))
#         g_temp = SimpleGraph(edges)
#         if isempty(edges)
#             add_vertices!(g, length(states))
#         else
#             g = union(g, g_temp)
#         end
#         i += length(states)

#         # Identify states up to identity generators
#         state_labels = vertex_labels(states, legend)
#         remove_vacuum!(state_labels, g)
#     end
#     push!(gg, g)
# end

# window_height = 10000
# window_width = 10000
# colors = ["paleturquoise", "chartreuse", "thistle1", "pink",
# "gold", "wheat", "olivedrab1", "palegreen", "turquoise1",
# "lightgreen", "plum1", "plum", "violet", "hotpink"]

# @pdf begin
#     background("grey10")
#     sethue("orange")
#     ng = length(gg)
#     N = convert(Int, ceil(sqrt(ng)))
#     tiles = Tiler(window_height, window_width, N, N)
#     setline(0.5)
#     for ((pos, n), (perm, _)) in zip(tiles, d)
#         @layer begin
#             n > ng && break
#             Karnak.translate(pos)
#             sethue("gray")
#             bbox = BoundingBox(box(O, tiles.tilewidth, tiles.tileheight, :stroke))
#             if nv(gg[n]) == 1
#                 layout = [Point(0.0,0.0)]
#             else
#                 layout = spring
#             end
#             sethue(colors[mod1(n, end)])
#             drawgraph(gg[n], layout = layout,
#             boundingbox = bbox,
#             vertexshapesizes = 4,
#             vertexfillcolors = colorant"red",
#             edgestrokeweights = 1)
#             sethue("cyan")
#             fontsize(14)
#             Karnak.text(string(Float64.(perm[1,:]))[2:end-1], halign = :center, boxbottomcenter(bbox) - (0, 20))
#             Karnak.text(string(Float64.(perm[2,:]))[2:end-1], halign = :center, boxbottomcenter(bbox) - (0, 5))
#         end
#     end
# end window_height+100 window_width+100 "pics/fragile/connectivity_full_merged_L$(L)_i$(interaction_range)"