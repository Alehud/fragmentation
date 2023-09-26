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


Plots.CURRENT_PLOT.nullableplot = nothing

# Spatial distribution of <n_b>
L = 180
n = 8
w = 3
data = deserialize("data/fragile/bs_L$(L)_n$(n)/bs_L$(L)_n$(n)_w$(w)_300.dat")
exp_vals = data["exp_vals_all"][1]
L = data["L"]
n = data["n"]
w = data["n_waves"]
t_therm = data["t_therm"]
t_meas = data["t_meas"]
plot([1:L;], exp_vals[[1, 75, 300]], lw=2, marker_z=[1:length(exp_vals);], palette=palette(:thermal, length(exp_vals[1:3])), colorbar=true, legend=false,
    colorbar_title=L"t/t_\mathrm{cycle}", colorbar_titlefontsize=14, labelfontsize=14)
title!(L"L=%$(L), n=%$(n), w=%$(w), t_\mathrm{cycle}=%$(int2label(t_meas))", titlefontsize=14)
xlabel!(L"x")
ylabel!(L"\left\langle n_b \right\rangle")
# savefig("pics/fragile/bs_expval_longword_L$(L)_n$(n)_w$(w)_$(int2str(t_meas))")


# <n_b>^2 over time
L = 180
n = 8
w = 3
data = deserialize("data/fragile/bs_L$(L)_n$(n)/bs_L$(L)_n$(n)_w$(w)_300.dat")
exp_vals = data["exp_vals_all"][1]
L = data["L"]
n = data["n"]
w = data["n_waves"]
t_therm = data["t_therm"]
t_meas = data["t_meas"]
tt = data["tt"][1] / t_meas
y_data = sum.((x -> x .^ 2).(exp_vals)) / L
plot([1:length(y_data);], y_data, legend=false, lw=2)
plot!([tt]; seriestype = :vline)
plot!([y_data[1]*0.15]; color=:lightgray, seriestype = :hline)
plot!([y_data[1]*0.10]; seriestype = :hline)
xlabel!(L"t")
ylabel!(L"$\left\langle n_b \right\rangle^2$ averaged over all sites")
title!(L"L=%$(L), n=%$(n), w=%$(w), t_\mathrm{cycle}=%$(int2label(t_meas))", titlefontsize=14)
# savefig("pics/fragile/bs_nb2_longword_L$(L)_n$(n)_w$(w)_$(int2str(t_meas))")


# Thermalization time (1/10) vs system size L (for fixed n=10)
means = Float64[]
errors = Float64[]
L_list = [45:5:130;]
for L in L_list
    data = deserialize("data/fragile/bs_n10/bs_L$(L)_n10_w1_10^4.dat")
    tt = data["tt"]
    tt = tt[tt.>0]
    println("L: $(L), experiments thermalized: $(length(tt))")
    push!(means, mean(tt))
    push!(errors, std(tt))
end
# linear_fit = Polynomials.fit(L_list, log.(means), 1)
# println(linear_fit.coeffs)
plot(L_list, means, yerror=errors,
    #  label=L"%$(round(exp(linear_fit.coeffs[1]); digits=2)) \cdot \exp(%$(round(linear_fit.coeffs[2]; digits=2))L)", 
    # yaxis=:log,
    # yticks=[10^4, 10^5, 10^6, 10^7, 10^8], 
    xticks=L_list,
    legend=false,
    legend_font=14,
    lw=2,
    labelfontsize=16)
xlabel!(L"$L$")
ylabel!(L"$t_{1/10}$")
title!(L"n=10", titlefontsize=16)
savefig("pics/fragile/bs_therm_time_n10_log")


# Thermalization time (1/10) vs system size L
means = Float64[]
errors = Float64[]
for L in 60:10:130
    data = deserialize("data/fragile/bs_L$(L)_n$(L÷10).dat")
    tt = data["tt"]
    push!(means, mean(tt))
    push!(errors, std(tt))
end
linear_fit = Polynomials.fit([60:10:130;], log.(means), 1)
println(linear_fit.coeffs)
plot([60:10:130;], means, yerror=errors,
    label=L"%$(round(exp(linear_fit.coeffs[1]); digits=2)) \cdot \exp(%$(round(linear_fit.coeffs[2]; digits=2))L)",
    yaxis=:log, yticks=[10^4, 10^5, 10^6, 10^7, 10^8], xticks=[60:10:130;],
    legend=:bottomright, legend_font=14, lw=2, labelfontsize=16)
xlabel!(L"$L$")
ylabel!(L"$t_{1/10}$")
title!(L"n=L/10", titlefontsize=16)
savefig("pics/fragile/bs_therm_time")



# Animation
L = 120
data = deserialize("data/fragile/bs_L$(L)_n$(L÷10)_one.dat")
exp_vals = data["exp_vals_all"][1]
nb2 = sum.((x -> x .^ 2).(exp_vals)) / L
n = data["n"]
t_therm = data["t_therm"]
t_meas = data["t_meas"]
anim = @animate for t in 1:100
    p1 = plot([1:L;], exp_vals[t], legend=false, lw=4, labelfontsize=32, xlabel=L"x", ylabel=L"\left\langle n_b \right\rangle", ylims=(-0.5, 0.55),
        title=L"t=%$(t),  t_\mathrm{cycle}=%$(int2label(t_meas))", titlefont=24, tickfontsize=20)
    p2 = plot([1:t;], nb2[1:t], legend=false, lw=4, labelfontsize=32, xlabel=L"t", ylabel=L"$\overline{\left\langle n_b \right\rangle^2}$",
        ylims=(0.0, 0.11), xlims=(0, 100), title=L"t=%$(t),  t_\mathrm{cycle}=%$(int2label(t_meas))", titlefont=24, tickfontsize=20)
    plot(p1, p2, windowsize=(2000, 800), margin=1.5 * Plots.cm)
end
gif(anim, "pics/fragile/anim_$(L).gif", fps=7)


# <nb>^2 over normalized time for different L, for a single experiment
Plots.CURRENT_PLOT.nullableplot = nothing
color_palette = palette(:rainbow, 8)
i = 1
for L in 60:10:130
    data = deserialize("data/fragile/bs_L$(L)_n$(L÷10)_one.dat")
    exp_vals = data["exp_vals_all"][1]
    t_meas = data["t_meas"]
    nb2 = sum.((x -> x .^ 2).(exp_vals)) / L
    tc = t_meas * findfirst(x -> x ≤ nb2[1] / 10, nb2)
    println(tc)
    plot!([1:t_meas:length(nb2)*t_meas;] / tc, nb2, legend=:topright, color=color_palette[i], lw=2, label=L"L=%$(L)", labelfontsize=13)
    i += 1
end
xlabel!(L"t/t_\mathrm{c}")
ylabel!(L"$\left\langle n_b \right\rangle^2$ averaged over all sites")
title!("single run")
xlims!((0.0, 2))
savefig("pics/fragile/bs_nb2_scaled_time_one")


# <nb>^2 over normalized time for different L, where <nb>^2 is averaged over many experiments
Plots.CURRENT_PLOT.nullableplot = nothing
color_palette = palette(:rainbow, 8)
i = 1
for L in 60:10:130
    data = deserialize("data/fragile/bs_L$(L)_n$(L÷10).dat")
    t_meas = data["t_meas"]
    experiments = data["experiments"]
    exp_vals_all = data["exp_vals_all"]
    nb2_all = (exp_vals -> sum.((x -> x .^ 2).(exp_vals)) / L).(exp_vals_all)
    tc_all = (nb2 -> t_meas * findfirst(x -> x ≤ nb2[1] / 10, nb2)).(nb2_all)
    t_all = [[1:length(nb2);] * t_meas / tc for (nb2, tc) in zip(nb2_all, tc_all)]
    lin_int_all = [linear_interpolation(t, nb) for (t, nb) in zip(t_all, nb2_all)]
    t = [t_meas/minimum(tc_all):0.001:1;]
    nb2_all_interpolated = [lin_int(t) for lin_int in lin_int_all]
    nb2_all_mean = mean(nb2_all_interpolated)
    plot!(t, nb2_all_mean, legend=:bottomleft, color=color_palette[i], lw=2, label=L"L=%$(L) \,\, (%$(experiments) \, \mathrm{runs})", labelfontsize=13)
    i += 1
end
xlabel!(L"t/t_\mathrm{c}")
ylabel!(L"$\left\langle n_b \right\rangle^2$ averaged over all sites")
title!("averaged over runs")
xlims!((0.0, 1.05))
savefig("pics/fragile/bs_nb2_scaled_time")