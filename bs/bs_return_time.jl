include("../src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Graphs
using Karnak
using Plots
using Colors
using LaTeXStrings
using StatsBase
using Polynomials
using Printf

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
L = 12
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


s_init = fill(Int8(0), L)
println("s_init: $s_init")

num_return_times = 10
t = @elapsed begin
    return_times = get_return_times(s_init, H, move!; num_returns=num_return_times)
end

println("$t sec")
returns_per_day = Int64(floor(24*60*60 / t))*num_return_times
println(returns_per_day, " estimated returns in a day")
flush(stdout)


return_times = get_return_times(s_init, H, move!; num_returns=returns_per_day)
serialize("bs/data/bs_return_times_L$(L)_$(length(return_times) ÷ 1000)k.dat", return_times)



L = 6
return_times = deserialize("bs/data/bs_return_times_L$(L).dat")
h = StatsBase.fit(Histogram, return_times, nbins=50)
edg = collect(h.edges[1])
x_data = (edg[1:end-1] .+ edg[2:end])/2
y_data = h.weights

histogram(return_times, bins=edg, label="Numerics ($(length(return_times) ÷ 1000)k returns)", yscale=:log10, legend=:topright, color=:gray, minorgrid=true)

linear_fit = Polynomials.fit(x_data, log.(replace(y_data, 0=>1)), 1)
a = @sprintf "%.2e" exp(linear_fit.coeffs[1])
b = @sprintf "%.2e" linear_fit.coeffs[2]

x13_fit = Polynomials.fit(x_data.^(1/3), log.(replace(y_data, 0=>1)), 1)
c = @sprintf "%.2f" x13_fit[1]
d = @sprintf "%.2e" x13_fit[2]

Plots.plot!(x_data, exp.(linear_fit.(x_data)), label=L"$ae^{bt}$", lw=3, color=:red, yscale=:log10, legend_font=10)
# Plots.plot!(x_data, exp(linear_fit.coeffs[1]) * exp.(linear_fit.coeffs[2]*x_data.^(1/3)), label=L"$ae^{bt^{1/3}}$", lw=3, color=:green, yscale=:log10)
title!("Return time probability distribution function \n Baumslag-Solitar group, L=$(L)", titlefontsize=10)
xlabel!(L"t")
ylabel!("# of returns")
savefig("pics/bs_return_times_L$(L)")