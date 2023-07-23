using IterTools
using Combinatorics

function N(n::Integer, label::Tuple{Vararg{Char}}, colors::Vector{Char}; pbc::Bool=false)
    N_list = Int8[]
    cols_list = NTuple{n, Char}[]
    for cols in product(fill(colors, n)...)
        if valid_label(cols, pbc=pbc)
            # println("cols: $cols")
            s = 0
            len = length(label)
            for i in product(fill(1:len, n)...)
                if all(i[1:end-1] .< i[2:end])
                    # println("i: $i")
                    if label[collect(i)] == cols
                        s += (-1)^sum(i)
                    end
                end
            end
            push!(N_list, s)
            push!(cols_list, cols)
        end
    end
    return N_list, cols_list
end

function N1(label::Union{Tuple{Vararg{Char}}, Vector{Char}}, colors::Vector{Char})
    N = Int8[]
    for α in colors
        # print("N $α: ")
        s = 0
        for i in eachindex(label)
            if label[i] == α
                s += (-1)^i
            end
        end
        # println(s)
        append!(N, s)
    end
    return N
end

function N2(label::Union{Tuple{Vararg{Char}}, Vector{Char}}, colors::Vector{Char})
    N = Int8[]
    for α in colors
        for β in colors
            if α ≠ β
                # print("N $α$β: ")
                s = 0
                for i in eachindex(label)
                    if label[i] == α
                        for j in i+1:length(label)
                            if label[j] == β
                                s += (-1)^(i+j)
                            end
                        end
                    end
                end
                # println(s)
                append!(N, s)
            end
        end
    end
    return N
end

function no_intersections(label::Union{Tuple{Vararg{Char}}, Vector{Char}})
    st = collect(label)
    while !isempty(st)
        flag = true
        for i in 1:length(st)-1
            if st[i] == st[i+1]
                deleteat!(st, i+1)
                deleteat!(st, i)
                flag = false
                break
            end
        end
        if flag
            break
        end
    end
    return isempty(st)
end

function valid_label(label::Union{Tuple{Vararg{Char}}, Vector{Char}}; pbc::Bool=false)
    st = collect(label)
    if pbc
        return all(circshift(st, 1) .≠ st)
    else
        return all(st[1:end-1] .≠ st[2:end])
    end
end

function unique_permutations(x::T, prefix=T()) where T
    if length(x) == 1
        return [[prefix; x]]
    else
        t = T[]
        for i in eachindex(x)
            if i > firstindex(x) && x[i] == x[i-1]
                continue
            end
            append!(t, unique_permutations([x[begin:i-1];x[i+1:end]], [prefix; x[i]]))
        end
        return t
    end
end


colors = ['R', 'G', 'B']
# colors = [0, 1]
# L = 16
# label = ('R', 'B', 'G', 'R', 'B', 'G', 'B', 'R', 'G', 'B', 'R', 'G')
# label = ('R', 'B', 'G', 'R', 'B', 'G', 'R', 'B', 'G', 'R', 'B', 'G', 'R', 'B', 'G', 'R', 'B', 'G')
# label = ('R', 'B', 'B', 'R', 'B', 'R', 'B', 'R', 'B', 'R', 'R', 'G')
# label = ('G', 'G', 'G', 'G', 'R', 'R')
# L = length(label)
# println("L: $L")
# P = N(1, label, colors)
# PP = N(2, label, colors)
# cols2 = PP[2]
# PPP = N(3, label, colors)
# println(P)
# println(PP)
# println(PPP)
# println(no_intersections(label))


## CHECKING THE RECURSION RELATION FOR 2-INDICES
# function N1_thru_N2(α, label, colors)
#     s = 0
#     pp = N(2, label, colors)
#     for β in colors
#         if β ≠ α
#             ind1 = findfirst(pp[2] .== [(α, β)])
#             ind2 = findfirst(pp[2] .== [(β, α)])
#             s += pp[1][ind1] - pp[1][ind2]
#         end
#     end
#     return s
# end

# function f(α, β)
#     s = 0
#     for j1 in 1:L÷2
#         s -= (label[2j1-1] == α)*(label[2j1] == β)
#         for γ in colors
#             if γ ≠ β
#                 pp = N(2, label[2j1+1:end], colors)
#                 ind1 = findfirst(pp[2] .== [(β, γ)])
#                 ind2 = findfirst(pp[2] .== [(γ, β)])
#                 s += (-(label[2j1-1] == α) + (label[2j1] == α))*(pp[1][ind1] - pp[1][ind2])
#             end
#         end
#     end
#     return s
# end

# # for (a, b) in cols2
# #     print(f(a, b), "  ")
# # end
# # println()

# println(label[L-1:end])
# println(N(2, label[L-1:end], colors))


# L = 10
# for label in product(fill(colors, L)...)
#     P1 = N2(label[[1:5; 10]], colors)
#     P2 = N2(label[5:10], colors)
#     n1 = count(x -> x ≠ 0, P1)
#     n2 = count(x -> x ≠ 0, P2)
#     if n1 == 0 && n2 == 0
#         PP = N2(label[[1:4; 6:9]], colors)
#         n = count(x -> x ≠ 0, PP)
#         if n ≠ 0
#             noint = no_intersections(label)
#             println(label, "  N2: ", PP, "  noint: ", noint)
#         end
#     end
# end


# # Check whether PP=0 implies P=0 (yes)
# for label in product(fill(colors, L)...)
#     P = N1(label, colors)
#     n1 = count(x -> x ≠ 0, P)
#     if n1 ≠ 0
#         PP = N2(label, colors)
#         n2 = count(x -> x ≠ 0, PP)
#         if n2 == 0
#             noint = no_intersections(label)
#             println(label, ":  N1: ", P, "  N2: ", PP, "  noint: ", noint)
#         end
#     end
# end


# # Check whether PPP=0 implies PP=0 (if P=0) (3 colors: yes up to L=12)
# for label in product(fill(colors, L)...)
#     P = N1(label, colors)
#     n1 = count(x -> x ≠ 0, P)
#     if n1 == 0
#         PP = N2(label, colors)
#         n2 = count(x -> x ≠ 0, PP)
#         if n2 ≠ 0
#             PPP, _ = N(3, label, colors)
#             n3 = count(x -> x ≠ 0, PPP)
#             if n3 == 0
#                 noint = no_intersections(label)
#                 println(label, ":  N1: ", N1(label, colors), "N2: ", PP, "  N3: ", PPP, "  noint: ", noint)
#             end
#         end
#     end
# end


# # Check whether P=0 and PP=0 implies no intersections
# for label in product(fill(colors, L)...)
#     # if 'R' in label && 'G' in label && 'B' in label && 'Y' in label && 'O' in label
#         P = N1(label, colors)
#         n1 = count(x -> x ≠ 0, P)
#         if n1 == 0
#             PP = N2(label, colors)
#             n2 = count(x -> x ≠ 0, PP)
#             if n2 == 0
#                 noint = no_intersections(label)
#                 if !noint
#                     println(label, ":  N1: ", P, "  N2: ", PP, "  noint: ", noint)
#                 end
#             end
#         end
#     # end
# end


for L in [4, 6, 8, 10, 12, 14, 16, 18]
    println("---------------------------------------------------------------------------------")
    println("L: $L")
    flag = false
    # Check whether P=0, PP=0 and PPP=0 implies no intersections (3 colors: yes up to L=16)
    for label in product(fill(colors, L)...)
        P = N1(label, colors)
        n1 = count(x -> x ≠ 0, P)
        if n1 == 0
            PP = N2(label, colors)
            n2 = count(x -> x ≠ 0, PP)
            if n2 == 0
                PPP, _ = N(3, tuple(label...), colors)
                n3 = count(x -> x ≠ 0, PPP)
                if n3 == 0
                    noint = no_intersections(label)
                    if !noint
                        flag = true
                        println(label, ":  N1: ", P, "  N2: ", PP, "  N3: ", PPP, "  noint: ", noint)
                    end
                end
            end
        end
    end
    if !flag
        println("For m=$(length(colors)) colors and L=$L, 1,2,3-index are enough for no intersection.")
    end
    flush(stdout)
end


println("END")