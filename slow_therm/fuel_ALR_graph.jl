include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using LinearAlgebra
using Graphs
using Luxor
using Karnak
using Colors
using Plots

function nice_string(s::Vector{Int8})
    return join(replace(s, 0 => '0', 1 => 'L', 2 => 'R', 3 => "B", 4 => 'A'))
end

function nice_print(s::Vector{Int8})
    println(nice_string(s))
end

function nice_print(states::Vector{Vector{Int8}})
    for s in states
        nice_print(s)
    end
end

# 0 = 0
# 1 = L
# 2 = R
# 3 = LR (B)
# 4 = A
dof_dim = 5
pbc = true

s_init = Int8[4,0,0,0,0]
println("s_init: $(s_init)")


# Construct Hamiltonian
L = length(s_init)
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
for site in 1:L-2+2*Int(pbc)
    # Creation of L fuel
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 4 && !(4 in s[2:3]) && !(isodd(s[2]) && isodd(s[3])), s -> Int8[4, (s[2]÷2)*2+1, (s[3]÷2)*2+1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 4 && isodd(s[2]) && isodd(s[3]), s -> Int8[4, s[2]-1, s[3]-1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 4 && isodd(s[2]) && isodd(s[3]), s -> Int8[4, s[2]-1, s[3]]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 4 && isodd(s[2]) && isodd(s[3]), s -> Int8[4, s[2], s[3]-1]))

    # Creation of R fuel
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[3] == 4 && !(4 in s[1:2]) && !(1 < s[1] < 4 && 1 < s[2] < 4), s -> Int8[(1-s[1]÷2)*2 + s[1], (1-s[2]÷2)*2 + s[2], 4]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[3] == 4 && 1 < s[1] < 4 && 1 < s[2] < 4, s -> Int8[s[1]-2, s[2]-2, 4]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[3] == 4 && 1 < s[1] < 4 && 1 < s[2] < 4, s -> Int8[s[1]-2, s[2], 4]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[3] == 4 && 1 < s[1] < 4 && 1 < s[2] < 4, s -> Int8[s[1], s[2]-2, 4]))

    # Movement of L fuel
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> isodd(s[1]) && isodd(s[2]) && (s[3] == 0 || s[3] == 2), s -> Int8[s[1]-1, s[2]-1, s[3]+1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> (s[1] == 0 || s[1] == 2) && (s[2] == 0 || s[2] == 2) && isodd(s[3]), s -> Int8[s[1]+1, s[2]+1, s[3]-1]))

    # Movement of R fuel
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] < 2 && 1 < s[2] < 4 && 1 < s[3] < 4, s -> Int8[s[1]+2, s[2]-2, s[3]-2]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> 1 < s[1] < 4 && s[2] < 2 && s[3] < 2, s -> Int8[s[1]-2, s[2]+2, s[3]+2]))

    # Movement of A
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> isodd(s[1]) && s[2] == 4 && 1 < s[3] < 4, s -> Int8[4, 0, s[3]]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> isodd(s[1]) && s[2] == 4 && 1 < s[3] < 4, s -> Int8[s[1], 0, 4]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 4 && s[2] == 0 && 1 < s[3] < 4, s -> Int8[1, 4, s[3]]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 4 && s[2] == 0 && 1 < s[3] < 4, s -> Int8[3, 4, s[3]]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> isodd(s[1]) && s[2] == 0 && s[3] == 4, s -> Int8[s[1], 4, 2]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> isodd(s[1]) && s[2] == 0 && s[3] == 4, s -> Int8[s[1], 4, 3]))
end
println(length(H_terms)÷2, " terms in the Hamiltonian")
H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);

# states_visited = Vector{Int8}[]
# creation_of_L = Tuple{Vector{Int8}, Vector{Int8}}[]
# creation_of_R = Tuple{Vector{Int8}, Vector{Int8}}[]
# movement_of_L = Tuple{Vector{Int8}, Vector{Int8}}[]
# movement_of_R = Tuple{Vector{Int8}, Vector{Int8}}[]
# movement_of_A = Tuple{Vector{Int8}, Vector{Int8}}[]
# for state in product(fill(Int8[0:dof_dim-1;], L)...)
#     s = collect(state)
#     push!(states_visited, s)
#     for (i, h) in enumerate(H_terms)
#         if h[3](s)
#             if !(h[4](s) in states_visited) 
#                 if 1 ≤ i ≤ 4
#                     push!(creation_of_L, (s, h[4](s)))
#                 elseif 5 ≤ i ≤ 8
#                     push!(creation_of_R, (s, h[4](s)))
#                 elseif 9 ≤ i ≤ 10
#                     push!(movement_of_L, (s, h[4](s)))
#                 elseif 11 ≤ i ≤ 12
#                     push!(movement_of_R, (s, h[4](s)))
#                 elseif 13 ≤ i ≤ 18
#                     push!(movement_of_A, (s, h[4](s)))
#                 else
#                     println("AAAAAAAAAAA")
#                 end
#             end
#         end
#     end
# end
# println("Creation of L fuel:")
# for tup in creation_of_L
#     println(nice_string(tup[1]), " <---> ", nice_string(tup[2]))
# end
# println("Creation of R fuel:")
# for tup in creation_of_R
#     println(nice_string(tup[1]), " <---> ", nice_string(tup[2]))
# end
# println("Movement of L fuel:")
# for tup in movement_of_L
#     println(nice_string(tup[1]), " <---> ", nice_string(tup[2]))
# end
# println("Movement of R fuel:")
# for tup in movement_of_R
#     println(nice_string(tup[1]), " <---> ", nice_string(tup[2]))
# end
# println("Movement of A particles:")
# for tup in movement_of_A
#     println(nice_string(tup[1]), " <---> ", nice_string(tup[2]))
# end


# Explore Krylov sector
states, ham = explore_connected_states(s_init, H, construct_ham=true, verbose=false);
println("$(length(states)) states found")
rows, cols, mels = findnz(ham)

# nice_print(states)
if length(states) == L*(dof_dim-1)^(L-1)
    println("Ergodic")
else
    println("Not ergodic. Missing states:")
    for s in product(fill(Int8[0:dof_dim-1;], L)...)
        if count(s .== 4) == 1
            if !(collect(s) in states)
                nice_print(collect(s))
            end
        end
    end
end


# Construct the graph
edges = Edge.(collect(zip(rows, cols)))
g = SimpleGraph(edges)
legend = [0 => "0", 1 => "L", 2 => "R", 3 => "B", 4 => "A"]
state_labels = vertex_labels(states, legend)
# remove_vacuum!(state_labels, g)
# println("$(length(state_labels)) states left after removing vacuum.")

# Plot the graph
y = [length(a_star(g, 1, v)) for v in vertices(g)]
H = 2000
W = 2000
# pts = vcat([between.(O + (-W/2, H/2 - d*H/maximum(y)), O + (W/2, H/2 - d*H/maximum(y)), range(0, 1, length=count(i -> i == d, y)+2))[2:end-1] 
#             for d in 0:maximum(y)]...)[invperm(sortperm(y))]
colors = [colorant"lightcoral", colorant"darkorange", colorant"gold", colorant"palegreen", colorant"skyblue", colorant"plum2", colorant"olivedrab", colorant"sandybrown", colorant"darkslateblue"]

@pdf begin
    background("grey10")
    sethue("orange")
    drawgraph(g, layout = stress,
    # vertexshapes = (vtx) -> box(O, 10(L+1), 20, :fill),
    # vertexlabels = state_labels,
    # vertexlabelfontsizes = 20,
    # vertexshapesizes = 12,
    # vertexlabeltextcolors = colorant"black",
    # vertexlabeloffsetdistances = 0,
    vertexshapesizes = 10,
    vertexfillcolors = (v) -> colors[findfirst('A', state_labels[v])],
    vertexstrokeweights = 0,
    edgelines = (k, s, d, f, t) ->
        Karnak.arrow(f, t, arrow_func(s, d, y),
            linewidth = 3,
            arrowheadlength = 0),
    edgestrokeweights = 4,
    edgecurvature = 30)
end H+100 W+100 "pics/slow_therm/fuel_ALR/L=$(L).pdf"
