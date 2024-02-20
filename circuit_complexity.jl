include("src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using Plots
using LaTeXStrings
using StatsBase


bit_lines = 3
gate_set_label = "Tof"

# Construct the gate set (a vector with functions)
# gate_set = Function[bits::Integer -> bits]    # Identity gate
# # NOT gates
# for i in 0:bit_lines-1
#     push!(gate_set, bits::Integer -> flip_bit(bits, i))
# end
# # Toffoli gates
# for i in 0:bit_lines-1, j in 0:bit_lines-1, k in j+1:bit_lines-1
#     if i ≠ j && i ≠ k
#         push!(gate_set, bits::Integer -> begin
#                                             if get_bit(bits, j) && get_bit(bits, k)
#                                                 return flip_bit(bits, i)
#                                             else
#                                                 return bits
#                                             end
#                                         end)
#     end
# end



# Calculate complexity of permutations and states that saturate the complexity
d_comp = Dict{Vector{Int8}, Integer}()                          # complexity of each permutation sector (length of the shortest circuit realizing this permutation)
d_comp_states = Dict{Vector{Int8}, Vector{Vector{Int8}}}()      # states that saturate the complexity (i.e., states with the number of non-identity gates equal to the complexity) for each permutation sector
identity_perm = Int8[0:(2^bit_lines - 1);]
d_comp[identity_perm] = 0
d_comp_states[identity_perm] = Vector{Int8}[[]]
Lmax = 7
for L in 1:Lmax
    println("L=$L -----------------------------------------------------------")
    flush(stdout)
    # d_noid = deserialize("data/circuit/d_noid_L$(L).dat")
    d_noid_len = deserialize("data/circuit_$(bit_lines)bits_$(gate_set_label)/d_noid_len/d_noid_len_L$(L).dat")
    
    # Update d_comp and d_comp_states with encountered group elements
    for g in keys(d_noid_len)
        if !(g in keys(d_comp))
            d_comp[g] = L
            # d_comp_states[g] = d_noid[g]
        end
    end
end
serialize("data/circuit_$(bit_lines)bits_$(gate_set_label)/d_comp.dat", d_comp)
# serialize("data/circuit/d_comp_states_until_L$Lmax.dat", d_comp_states)



# Calculate ω(K,N) - the number of N-gate circuits (without identities) with complexity K.
# ω[k+1, n] is ω(k, n) (note offset by 1, because complexity can be 0)
d_comp = deserialize("data/circuit/d_comp.dat")
N = 200
ω = zeros(BigInt, maximum(values(d_comp))+1, N)
d_noid_len = deserialize("data/circuit/d_noid_len_until_L200.dat")
for n in 1:N
    println("n: ", n)
    flush(stdout)
    for (p, num) in d_noid_len
        ω[d_comp[p]+1, n] += num[n]
    end
end
serialize("data/circuit/omega_until_N$N.dat", ω)


# Calculate ν(K) - the number of permutations with complexity K.
d_comp = deserialize("data/circuit_3bits_Tof/d_comp.dat")
nu = [count(x -> x == i, values(d_comp)) for i in 0:12]
serialize("data/circuit_3bits_Tof/nu.dat", nu)


# ω_list(K,N). Same matrix as ω(K,N), but every element is a vector with numbers of length-N circuits with complexity K in different permutation sectors 
d_comp = deserialize("data/circuit/d_comp.dat")
N = 200
ω_list = [Vector{BigInt}() for _ in 1:maximum(values(d_comp))+1, _ in 1:N]
d_noid_len = deserialize("data/circuit/d_noid_len_until_L200.dat")
for n in 1:N
    println("n: ", n)
    flush(stdout)
    for (p, num) in d_noid_len
        push!(ω_list[d_comp[p]+1, n], num[n])
    end
end
serialize("data/circuit/omega_list_until_N$(N).dat", ω_list)



N = 12
ν = deserialize("data/circuit/nu.dat")
ω_list = deserialize("data/circuit/omega_until_N200.dat")
log_nu = log2.(ν)
S_bar_1 = Real[]
S_bar_2 = Real[]
for K in 0:N
    sorted = sort(ω_list[(K,N)])
    push!(S_bar_1, log2(mean(sorted)))
    push!(S_bar_2, log2(maximum(sorted)))
end

x_nu = [0:9;]
x_S_bar = [0:1:12;]
nu_fit = Polynomials.fit(x_nu, log_nu[x_nu .+ 1], 1)
S_bar_fit_1 = Polynomials.fit(x_S_bar, S_bar_1[x_S_bar .+ 1], 1)
S_bar_fit_2 = Polynomials.fit(x_S_bar, S_bar_2[x_S_bar .+ 1], 1)

x = [0:12;]
plot(x, log_nu, seriestype=:scatter, label=L"$\log_2(\nu)$", color=:red, windowsize=(1200, 800), legend_font=14)
plot!(x, S_bar_1, seriestype=:scatter, label=L"$\bar{S}_1$", color=:blue)
plot!(x, S_bar_2, seriestype=:scatter, label=L"$\bar{S}_2$", color=:green)
plot!(x, nu_fit.(x), label="$(nu_fit)", color=:red, lw=2)
plot!(x, S_bar_fit_1.(x), label="$(S_bar_fit_1)", color=:blue, lw=2)
plot!(x, S_bar_fit_2.(x), label="$(S_bar_fit_2)", color=:green, lw=2)
xlabel!(L"Κ")
# # savefig("pics/circuit/alpha_beta")









L = 12
Ω = deserialize("data/circuit_$(bit_lines)bits_$(gate_set_label)/d_noid_len/d_noid_len_L$(L).dat")
comp = deserialize("data/circuit_$(bit_lines)bits_$(gate_set_label)/d_comp.dat")
Ω_comp = Dict([p => (Ω[p], comp[p]) for p in keys(Ω)])
x_data = [x[2] for x in collect(values(Ω_comp))]
y_data = [x[1] for x in collect(values(Ω_comp))]
if gate_set_label == "NotTof"
    title = "$(bit_lines) bit lines; Gate set: NOT + Toffoli (one polarity)"
elseif gate_set_label == "Tof"
    title = "$(bit_lines) bit lines; Gate set: bipolar Toffolis"
end
plot(x_data, log2.(y_data),
     seriestype=:scatter,
     markersize=6,
     markercolor=:red,
     labelfontsize=18,
     tickfontsize=14,
     titlefontsize=16,
     legendfontsize=14,
     title=title,
     label=L"L=%$(L)",
     legend=:bottomleft,
     ylim=(0,maximum(log2.(y_data))+1),
     xticks=[0:2:12;],
     windowsize=(1000, 600),
     margin=0.4 * Plots.cm
     )
xlabel!(L"Κ")
ylabel!(L"\log_2(\Omega)")
savefig("pics/circuit/S_vs_K_$(gate_set_label)")




nu = deserialize("data/circuit_$(bit_lines)bits_$(gate_set_label)/nu.dat")
if gate_set_label == "NotTof"
    title = "$(bit_lines) bit lines; Gate set: NOT + Toffoli (one polarity)"
elseif gate_set_label == "Tof"
    title = "$(bit_lines) bit lines; Gate set: bipolar Toffolis"
end
plot([0:length(nu)-1;], log2.(nu),
     seriestype=:scatter,
     markersize=6,
     markercolor=:blue,
     labelfontsize=18,
     tickfontsize=14,
     titlefontsize=16,
     legendfontsize=16,
     title=title,
     label=false,
     legend=:topleft,
    #  ylim=(0,maximum(log2.(Ω))+1),
     xticks=[0:12;],
     windowsize=(1000, 600),
     margin=0.4 * Plots.cm
     )
xlabel!(L"K")
ylabel!(L"\log_2(\nu)")
savefig("pics/circuit/nu_vs_K_$(gate_set_label)")



d_noid_len = deserialize("data\\circuit_$(bit_lines)bits_$(gate_set_label)\\d_noid_len\\d_noid_len_until_L100.dat")
if gate_set_label == "NotTof"
    title = "$(bit_lines) bit lines; Gate set: NOT + Toffoli (one polarity)"
elseif gate_set_label == "Tof"
    title = "$(bit_lines) bit lines; Gate set: bipolar Toffolis"
end
L_list = [2:18;]
plot([factorial(2^bit_lines)], 
      seriestype=:hline, 
      c=:black,
      lw=2,
      linestyle=:dash,
      label=L"(2^3)!",
      labelfontsize=18,
      tickfontsize=14,
      titlefontsize=16,
      xlabel=L"s",
      ylabel="# of permutations",
      title=title,
      legend=:outerright,
      legend_font=12,
      xticks=[10^n for n in 0:15],
      yticks=[10^n for n in 0:5],
      xaxis=:log,
      yaxis=:log,
      palette=palette(:darktest, length(L_list)+2),
      windowsize=(1000, 600),
      margin=0.4 * Plots.cm
)
for L in L_list
    num_perms = [x[L] for x in values(d_noid_len)]
    num_perms = num_perms[num_perms .> 0]
    cdf = ecdf(num_perms)

    plot!(x -> length(num_perms)*(1-cdf(x-1)), minimum(num_perms), maximum(num_perms)-1,
        label=L"L=%$(L)",
        lw=2,
        )
end
plot!()
savefig("pics/circuit/num_perms_vs_s_$(gate_set_label)")

