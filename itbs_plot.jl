include("src/Fragmentation.jl")
using .Fragmentation
using SparseArrays
using Serialization
using Plots
using Colors
using LaTeXStrings
using Statistics
using Polynomials
using Interpolations
using RollingFunctions
using StatsBase


Plots.CURRENT_PLOT.nullableplot = nothing

# Spatial distribution of <n_c>
L = 100
n = 4
data = deserialize("data/itbs/itbs_n$(n)/itbs_L$(L)_n$(n)_10^3.dat")
e = 1
exp_vals = data["exp_vals_all"][e]
L = data["L"]
n = data["n"]
w = data["n_waves"]
t_therm = data["t_therm"]
t_meas = data["t_meas"]
ind_to_show = 1:length(exp_vals)
plot([1:L;], exp_vals[ind_to_show], lw=2, marker_z=[1:length(exp_vals[ind_to_show]);], palette=palette(:thermal, length(exp_vals[ind_to_show])), colorbar=true, legend=false,
    colorbar_title=L"t/t_\mathrm{cycle}", colorbar_titlefontsize=14, labelfontsize=14)
title!(L"L=%$(L), n=%$(n), w=%$(w), t_\mathrm{cycle}=%$(int2label(t_meas))", titlefontsize=14)
xlabel!(L"x")
ylabel!(L"\left\langle n_c \right\rangle")
# savefig("pics/bs/bs_expval_longword_L$(L)_n$(n)_w$(w)_$(int2str(t_meas))")


# <n_c>^2 over time
L = 90
n = 4
data = deserialize("data/itbs/itbs_n$(n)/itbs_L$(L)_n$(n)_100.dat")
e = 1
exp_vals = data["exp_vals_all"][e]
L = data["L"]
n = data["n"]
w = data["n_waves"]
t_therm = data["t_therm"]
t_meas = data["t_meas"]
tt = data["tt"][e] / t_meas
y_data = sum.((x -> x .^ 2).(exp_vals)) / L
roll_window = 30
smooth_y_data = rollmean(y_data, roll_window)
plot([1:length(y_data);], y_data, legend=false, lw=2)
# plot!([tt]; seriestype = :vline)
# plot!([0.08]; color=:orange, seriestype = :hline)
# plot!([y_data[1]*0.15]; color=:lightgoldenrod, seriestype = :hline)
plot!([y_data[1]*0.10]; color=:lightgreen, seriestype = :hline)
plot!([1:length(smooth_y_data);] .+ (roll_window÷2), smooth_y_data, legend=false, lw=2, color="red")
xlabel!(L"t")
ylabel!(L"$\left\langle n_c \right\rangle^2$ averaged over all sites")
title!(L"L=%$(L), n=%$(n), w=%$(w), t_\mathrm{cycle}=%$(int2label(t_meas))", titlefontsize=14)
savefig("pics/itbs/itbs_nc2_L$(L)_n$(n)_$(int2str(t_meas))")


# Thermalization time (1/10) vs system size L (for fixed n)
means = Float64[]
errors = Float64[]
n = 4
L_list = [50:5:65; 66:70; 75:5:150;]
for L in L_list
    data = deserialize("data/itbs/itbs_n$(n)/only_tt/itbs_L$(L)_n$(n)_100.dat")
    tt = data["tt"]
    tt = tt[tt.>0]
    println("L: $(L), experiments thermalized: $(length(tt)) out of $(length(data["tt"]))")
    push!(means, mean(tt))
    push!(errors, std(tt))
end
# linear_fit = Polynomials.fit(L_list, log.(means), 1)
# println(linear_fit.coeffs)
plot(L_list, means / 10^4, yerror=errors / 10^4,
    seriestype=:scatter,
    markersize=6,
    color=:orange,
    markerstrokewidth=3,
    #  label=L"%$(round(exp(linear_fit.coeffs[1]); digits=2)) \cdot \exp(%$(round(linear_fit.coeffs[2]; digits=2))L)", 
    # yaxis=:log,
    # yticks=[10^3, 5*10^3, 10^4, 10^5], 
    xticks=[20:20:150;],
    legend=false,
    labelfontsize=24,
    tickfontsize=18,
    margin=0.3 * Plots.cm
    )
xlabel!(L"$L$")
ylabel!(L"$t_\mathrm{th} / 10^4$")
# title!(L"n=%$(n)", titlefontsize=16)
savefig("pics/itbs/itbs_t_th_n$(n).pdf")


# Expansion length vs n fow w_large(n)
n_list = [3:5;]
EL = [37,69,230]
linear_fit = Polynomials.fit(n_list, log2.(EL), 1)
println(linear_fit.coeffs)
plot(n_list, 2.17749*(2).^(1.31 * n_list),
      yaxis=:log2,
      label=L"2.18 \cdot 2^{1.31 n}",
      color=:blue,
      legend=:topleft,
      legend_font=26,
      lw=4, 
      labelfontsize=26,
      tickfontsize=20)
plot!(n_list, EL,
    seriestype=:scatter,
    markersize=8,
    markerstrokewidth=4,
    color=:red,
    yticks=[32, 64, 128, 256, 512], 
    xticks=n_list,
    label=false)
xlabel!(L"$n$")
ylabel!(L"$\mathrm{EL}$")
# title!(L"n=L/10", titlefontsize=22)
savefig("pics/itbs/itbs_expansion_length.pdf")


# Animation
L = 90
n = 4
data = deserialize("data/itbs/itbs_n$(n)/itbs_L$(L)_n$(n)_100.dat")
exp_vals = data["exp_vals_all"][1]
nb2 = sum.((x -> x .^ 2).(exp_vals)) / L
n = data["n"]
t_therm = data["t_therm"]
t_meas = data["t_meas"]
t_max = length(nb2)
anim = @animate for t in 1:20:t_max+600
    p1 = plot([1:L;], exp_vals[min(t, t_max)], legend=false, lw=4, labelfontsize=36, xlabel=L"x", ylabel=L"\left\langle n_\mathtt{c} \right\rangle", ylims=(-0.7, 0.65),
        title=L"L=%$(L), n=%$(n),  T=%$(int2label(t_meas)), t/T=%$(t)", titlefont=30, tickfontsize=24)
    p2 = plot([1:min(t, t_max);], nb2[1:min(t,t_max)], legend=false, lw=4, labelfontsize=36, xlabel=L"t/T", ylabel=L"$\overline{\left\langle n_\mathtt{c} \right\rangle^2}$",
        ylims=(0.0, maximum(nb2)*1.02), xlims=(0, t_max), title=L"L=%$(L), n=%$(n),  T=%$(int2label(t_meas)), t/T=%$(t)", titlefont=30, tickfontsize=24)
    plot(p1, p2, windowsize=(2000, 800), margin=1.5 * Plots.cm)
end
gif(anim, "pics/itbs/anim_$(L).gif", fps=10)