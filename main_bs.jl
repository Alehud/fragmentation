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

e = Int64[1 0; 0 1]
a = Int64[1 1; 0 1]
ā = Int64[1 -1; 0 1]
b = Int64[2 0; 0 1]
b̄ = [1//2 0; 0 1]

interaction_range = 3

# Construct the generator set (a vector with matrices)
generator_set = Matrix{Union{Int64, Rational{Int64}}}[e, a, b, ā, b̄]
dof_dim = length(generator_set)


# Construct a dictionary, where keys are computations (permutations) and values are vectors of states realizing this computation
d = Dict{Matrix{Union{Int64, Rational{Int64}}}, Vector{Vector{Int8}}}()
d_len = Dict{Matrix{Union{Int64, Rational{Int64}}}, Integer}()
identity_states = Vector{Int8}[]
for state in product(fill([0:(dof_dim-1);], interaction_range)...)
    group_elem = copy(e)
    for s in collect(state)
        group_elem *= generator_set[s+1]
    end
    d[group_elem] = push!(get(d, group_elem, Vector{Int8}[]), collect(state))
    d_len[group_elem] = get(d_len, group_elem, 0) + 1
    if group_elem == e
        push!(identity_states, collect(state))
    end
end
# findfirst(states -> [0,0,0] in states, d)     # how to find a permutation realized by a specific state
# filter!(p -> length(p.second) > 1, d)   # Remove states that are not flippable
# flippable_states = collect(values(d))    # leave only values (states)

# serialize("data/fragile/d_i$(interaction_range).dat", d)
# serialize("data/fragile/d_len_i$(interaction_range).dat", d_len)
# serialize("data/fragile/identity_states_i$(interaction_range).dat", identity_states)


# Construct Hamiltonian
L = 7
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
H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);


d = deserialize("data/fragile/d_i$(L).dat")
gg = SimpleGraph{Int64}[]
legend = ["0" => "_", "1" => "a", "2" => "b", "3" => "ā", "4" => "b̄"]
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

        # Identify states up to identity generators
        state_labels = vertex_labels(states, legend)
        remove_vacuum!(state_labels, g)
    end
    push!(gg, g)
end



window_height = 10000
window_width = 10000

# y = [length(a_star(g, 1, v)) for v in vertices(g)]
# pts = vcat([between.(O + (-window_width/2, window_height/2 - d*window_height/maximum(y)), O + (window_width/2, window_height/2 - d*window_height/maximum(y)), range(0, 1, length=count(i -> i == d, y)+2))[2:end-1] 
#             for d in 0:maximum(y)]...)[invperm(sortperm(y))]
# Stress(initialpos = [(pt.x, pt.y) for pt in pts], iterations = 100)


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
            Karnak.text(string(Float64.(perm[1,:]))[2:end-1], halign = :center, boxbottomcenter(bbox) - (0, 20))
            Karnak.text(string(Float64.(perm[2,:]))[2:end-1], halign = :center, boxbottomcenter(bbox) - (0, 5))
        end
    end
end window_height+100 window_width+100 "pics/fragile/connectivity_full_merged_L$(L)_i$(interaction_range)"