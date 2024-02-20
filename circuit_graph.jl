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

bit_lines = 3
interaction_range = 3
N_computations = factorial(2^bit_lines)
dof_dim = 7

d = deserialize("data/circuit/d_L$(interaction_range).dat")


# Construct Hamiltonian
L = 5
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
for site in 1:(L-(interaction_range-1))
    for (computation, states) in d
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


# Explore Hilbert space from some initial state
s_init = fill(Int8(0), L)
println("s_init: $s_init")
states, ham = explore_connected_states(s_init, H, construct_ham=true)
println("$(length(states)) states found")
# println("states:")
# for state in states
#     println(state)
# end
# println("ham: ", ham)
display(ham)
rows, cols, mels = findnz(ham);
# serialize("data/circuit/states_connected_to_identity_i$(interaction_range)_L$(L).dat", states)

# Construct the graph
edges = Edge.(collect(zip(rows, cols)))
g = SimpleGraph(edges)
legend = ["0" => "_", "1" => "1", "2" => "2", "3" => "3", "4" => "4", "5" => "5", "6" => "6"]
state_labels = vertex_labels(states, legend)
remove_vacuum!(state_labels, g)
println("$(length(state_labels)) states left after removing vacuum.")

# Plot the graph
y = [length(a_star(g, 1, v)) for v in vertices(g)]
H = 2000
W = 2000
pts = vcat([between.(O + (-W/2, H/2 - d*H/maximum(y)), O + (W/2, H/2 - d*H/maximum(y)), range(0, 1, length=count(i -> i == d, y)+2))[2:end-1] 
            for d in 0:maximum(y)]...)[invperm(sortperm(y))]
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
        Karnak.arrow(f, t, arrow_func(s, d, y),
            linewidth = 3,
            arrowheadlength = 0),
    edgestrokeweights = 4,
    edgecurvature = 30)
end H+100 W+100 "pics/circuit/circuit_$(L)x$(bit_lines)_i$(interaction_range)"


# Construct the full connectivity graph for the whole Hilbert space
d = deserialize("data/circuit/d_L$(L).dat")
gg = SimpleGraph{Int64}[]
for (perm, states_symm_sector) in d
    # Connectivity within symmetry sector
    states_symm_sector_list, hams = explore_full_space(H, states_symm_sector; construct_ham=true)
    g = SimpleGraph()
    i = 0
    for (states, ham) in zip(states_symm_sector_list, hams)
        rows, cols, mels = findnz(ham);
        rows .+= i
        cols .+= i
        edges = Edge.(collect(zip(rows, cols)))
        g_temp = SimpleGraph(edges)
        if isempty(edges)
            add_vertices!(g, length(states))
        else
            g = union(g, g_temp)
        end
        i += length(states)
    end
    push!(gg, g)
end
d_comp = deserialize("data/circuit/d_comp.dat")


# Plot the graph (separated into boxes) of the full Hilbert space
colors = ["paleturquoise", "chartreuse", "thistle1", "pink",
"gold", "wheat", "olivedrab1", "palegreen", "turquoise1",
"lightgreen", "plum1", "plum", "violet", "hotpink"]
window_height = 10000
window_width = 10000
@pdf begin
    background("grey10")
    sethue("orange")
    ng = length(gg)
    N = convert(Int, ceil(sqrt(ng)))
    tiles = Tiler(window_height, window_width, N, N)
    setline(0.5)
    for ((pos, n), (perm, _)) in zip(tiles, d)
        @layer begin
            n > ng && break
            Karnak.translate(pos)
            sethue("gray")
            bbox = BoundingBox(box(O, tiles.tilewidth, tiles.tileheight, :stroke))
            if nv(gg[n]) == 1
                layout = [Point(0.0,0.0)]
            else
                layout = spring
            end
            sethue(colors[mod1(n, end)])
            drawgraph(gg[n], layout = layout,
            boundingbox = bbox,
            vertexshapesizes = 4,
            vertexfillcolors = colorant"red",
            edgestrokeweights = 1)
            sethue("cyan")
            fontsize(14)
            Karnak.text(string(perm, "           Κ = ", d_comp[perm]), halign = :center, boxbottomcenter(bbox)-5)
        end
    end
end window_height+100 window_width+100 "pics/circuit/circuit_full_$(L)x$(bit_lines)_i$(interaction_range)"
