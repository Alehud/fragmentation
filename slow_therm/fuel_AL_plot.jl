include("../src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Plots
using Colors
using LaTeXStrings
using Statistics
using Polynomials
using RollingFunctions
using StatsBase



Plots.CURRENT_PLOT.nullableplot = nothing


L = 12
T = 10
data = deserialize("data/slow_therm/fuel_AL/L$(L)_T$(T)_long.dat")
evals_all = data["evals_all"]
exper = 3
evals = evals_all[exper]
plot([1:length(evals);], evals, 
     legend=false, 
     lw=3,
     labelfontsize=22,
     tickfontsize=12,
     margin=0.3 * Plots.cm)
xlabel!(L"t/T")
ylabel!("Position of A")
title!(L"L=%$(L), T=%$(int2label(T))", titlefontsize=20)
# savefig("pics/slow_therm/.pdf")



means = Float64[]
errors = Float64[]
L_list = [10:20;]
T_list = [3,5,10,30,50,100,300,500,1000,3000,5000]
for (L, T) in zip(L_list, T_list)
     data = deserialize("data/slow_therm/fuel_AL/L$(L)_T$(T).dat")
     evals_all = data["evals_all"]
     t_cr = (y -> findfirst(x -> x ≠ L÷2+1, y)).(evals_all)
     t_cr = t_cr[.!(isnothing.(t_cr))]*T
     println("L: $(L), experiments thermalized: $(count(x -> 0 < x, t_cr)), mean: $(mean(t_cr) / T)*T")
     push!(means, mean(t_cr))
     push!(errors, std(t_cr))
end
plot(L_list, means, yerror=errors,
     seriestype=:scatter,
     yaxis=:log,
     markersize=6,
     markerstrokewidth=2,
     color=:violet,
     yticks=[10^3, 10^4, 10^5, 10^6, 10^7, 10^8], 
     xticks=L_list,
     legend_font=14,
     lw=2,
     label=false,
     labelfontsize=22,
     tickfontsize=16,
     margin=0.3 * Plots.cm
)
linear_fit = Polynomials.fit(L_list, log2.(means), 1)
println(linear_fit.coeffs)
plot!(L_list, 2^(linear_fit.coeffs[1]) * 2 .^ (linear_fit.coeffs[2]*L_list),
      lw=2,
      c=:blue,
      label=L"0.458 \cdot 2^{1.209 L}",
      legend=:bottomright,
      )
xlabel!(L"L")
ylabel!(L"t_\mathrm{cr}")
# title!(L"", titlefontsize=20)
savefig("pics/slow_therm/fuel_AL/averaged/tcr(L).pdf")

