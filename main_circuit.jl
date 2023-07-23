include("src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using Graphs
using Karnak
using Colors


bit_lines = 3
interaction_range = 3
N_computations = factorial(2^bit_lines)


# Construct the gate set (a vector with functions)
gate_set = Any[bits::Integer -> bits]    # Identity gate
# NOT gates
for i in 0:bit_lines-1
    push!(gate_set, bits::Integer -> flip_bit(bits, i))
end
# Toffoli gates
for i in 0:bit_lines-1, j in 0:bit_lines-1, k in j+1:bit_lines-1
    if i ≠ j && i ≠ k
        push!(gate_set, bits::Integer -> begin
                                            if get_bit(bits, j) && get_bit(bits, k)
                                                return flip_bit(bits, i)
                                            else
                                                return bits
                                            end
                                        end)
    end
end
dof_dim = length(gate_set)


# Construct a dictionary, where keys are computations (permutations) and values are vectors of states realizing this computation
d = Dict{Vector{<:Integer}, Vector{Vector{<:Integer}}}()
d_len = Dict{Vector{<:Integer}, Integer}()
vacuum_states = Vector{<:Integer}[]
for state in [collect(x) for x in product(fill([0:(dof_dim-1);], interaction_range)...)]
    bits = [0:2^bit_lines-1;]
    for s in state
        bits = map(gate_set[s+1], bits)
    end
    d[bits] = push!(get(d, bits, Vector{<:Integer}[]), state)
    d_len[bits] = get(d_len, bits, 0) + 1
    if bits == [0:2^bit_lines-1;]
        push!(vacuum_states, state)
    end
end
# findfirst(states -> [0,0,0] in states, d)     # how to find a permutation realized by a specific state
filter!(p -> length(p.second) > 1, d)   # Remove states that are not flippable
flippable_states = collect(values(d))    # leave only values (states)

# open("data/circuit/vacuum_states_i$(interaction_range).txt", "w") do file
#     println(file, "$(length(d[[0:7;]])) states connected to the vacuum") 
#     for state in d[[0:7;]]
#     println(file, state)
#     end
# end

# serialize("data/circuit/dict_i$(interaction_range).dat", d)
# serialize("data/circuit/dict_len_i$(interaction_range).dat", d_len)
# serialize("data/circuit/vacuum_states_i$(interaction_range).dat", vacuum_states)


L = 5
H_terms = Tuple[]
for site in 1:(L-(interaction_range-1))
    for (computation, states) in d
        # println(computation)
        # println("     ", states)
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

# s_init = fill(Int8(0), L)
# println("s_init: $s_init")
# states, ham = explore_connected_states(s_init, H, construct_ham=true)
# println("$(length(states)) states found")
# # println("states:")
# # for state in states
# #     println(state)
# # end
# # println("ham: ", ham)
# display(ham)
# rows, cols, mels = findnz(ham);
# serialize("data/circuit/states_connected_to_vacuum_i$(interaction_range)_L$(L).dat", states)



d_i = deserialize("data/circuit/dict_i$(L).dat")
gg = SimpleGraph{Int64}[]
for (perm, states_symm_sector) in d_i
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



window_height = 10000
window_width = 10000
# y = [length(a_star(g, 1, v)) for v in vertices(g)]
# pts = vcat([between.(O + (-window_width/2, window_height/2 - d*window_height/maximum(y)), O + (window_width/2, window_height/2 - d*window_height/maximum(y)), range(0, 1, length=count(i -> i == d, y)+2))[2:end-1] 
#             for d in 0:maximum(y)]...)[invperm(sortperm(y))]
# # Stress(initialpos = [(pt.x, pt.y) for pt in pts], iterations = 100)


# @pdf begin
#     background("grey10")
#     sethue("orange")
#     drawgraph(g, layout = spring,
#     vertexshapes = (vtx) -> box(O, 10(L+1), 20, :fill),
#     vertexlabels = state_labels,
#     vertexlabelfontsizes = 20,
#     vertexshapesizes = 12,
#     vertexlabeltextcolors = colorant"black",
#     vertexlabeloffsetdistances = 0,
#     edgelines = (k, s, d, f, t) ->
#         arrow(f, t, arrow_func(s, d, y),
#             linewidth = 3,
#             arrowheadlength = 0),
#     edgestrokeweights = 4,
#     edgecurvature = 30)
# end window_height+100 window_width+100 "pics/circuit/circuit_$(L)x$(bit_lines)_i$(interaction_range)"

colors = ["paleturquoise", "chartreuse", "thistle1", "pink",
"gold", "wheat", "olivedrab1", "palegreen", "turquoise1",
"lightgreen", "plum1", "plum", "violet", "hotpink"]

@pdf begin
    background("grey10")
    sethue("orange")
    ng = length(gg)
    N = convert(Int, ceil(sqrt(ng)))
    tiles = Tiler(window_height, window_width, N, N)
    setline(0.5)
    for ((pos, n), (perm, _)) in zip(tiles, d_i)
        @layer begin
            n > ng && break
            translate(pos)
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
            text(string(perm), halign = :center, boxbottomcenter(bbox)-5)
        end
    end
end window_height+100 window_width+100 "pics/circuit/circuit_full_$(L)x$(bit_lines)_i$(interaction_range)"

# @pdf begin
#     background("grey10")
#     sethue("orange")
#     drawgraph(gg[165], layout = [Point(0.0, 0.0), Point(1000.0, 0.0)],
#     vertexshapes = :square,
#     vertexshapesizes = 40,
#     vertexfillcolors = colorant"red",
#     edgestrokeweights = 2)
# end window_height+100 window_width+100 "pics/circuit/circuit_$(L)x$(bit_lines)_i$(interaction_range)"


