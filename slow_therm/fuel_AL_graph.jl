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

function nice_print(s::Vector{Int8})
    println(join(replace(s, 0 => '0', 1 => 'L', 2 => 'A')))
end

function nice_print(states::Vector{Vector{Int8}})
    for s in states
        nice_print(s)
    end
end

# 0 = 0
# 1 = L
# 2 = A
dof_dim = 3
pbc = true

s_init = Int8[2,0,0,0,0,0,0,0,0]
println("s_init: $(s_init)")


# Construct Hamiltonian
L = length(s_init)
H_terms = Tuple{Real, Vector{<:Integer}, Function, Function}[]
for site in 1:L-2+2*Int(pbc)
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 1 && s[2] == 2 && s[3] ≠ 2, s -> Int8[1,0,2]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[1,0,2], s -> Int8[1,2,0]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[1,0,2], s -> Int8[1,2,1]))

    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s[1] == 2 && s[2] ≠ 2 && s[3] ≠ 2 && s ≠ Int8[2,1,1], s -> Int8[2,1,1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[2,1,1], s -> Int8[2,0,0]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[2,1,1], s -> Int8[2,0,1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[2,1,1], s -> Int8[2,1,0]))

    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[1,1,0], s -> Int8[0,0,1]))
    push!(H_terms, (1.0, [site, mod((site+1)-1, L)+1, mod((site+2)-1, L)+1], s -> s == Int8[0,0,1], s -> Int8[1,1,0]))
end
println(length(H_terms)÷2, " terms in the Hamiltonian")
H = Hamiltonian(dof_dim, H_terms, check_hermitian=true);


# Explore Krylov sector
states, ham = explore_connected_states(s_init, H, construct_ham=true, verbose=false);
println("$(length(states)) states found")
rows, cols, mels = findnz(ham)

# nice_print(states)
if length(states) == L*2^(L-1)
    println("Ergodic")
else
    println("Not ergodic. Missing states:")
    for s in product(fill(Int8[0,1,2], L)...)
        if count(s .== 2) == 1
            if !(collect(s) in states)
                nice_print(collect(s))
            end
        end
    end
end


# Construct the graph
edges = Edge.(collect(zip(rows, cols)))
g = SimpleGraph(edges)
legend = [0 => "0", 1 => "L", 2 => "A"]
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
end H+100 W+100 "pics/slow_therm/fuel_AL/L=$(L).pdf"
