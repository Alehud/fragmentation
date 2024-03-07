include("../src/Fragmentation.jl")
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

# Spatial distribution of <n_b>
L = 120
n = 12
w = 1
data = deserialize("data/bs/bs_one_runs/bs_L$(L)_n$(n)_one.dat")
e = 1
exp_vals = data["exp_vals_all"][e]
L = data["L"]
n = data["n"]
# w = data["n_waves"]
t_therm = data["t_therm"]
t_meas = data["t_meas"]
t_therm = data["tt"]
ind_to_show = 1:length(exp_vals)
plot([1:L;], exp_vals[ind_to_show], 
     lw=2, 
     marker_z=[0:length(exp_vals[ind_to_show]);], 
     palette=palette(:thermal, length(exp_vals[ind_to_show])), 
     legend=false,
     colorbar=true,
     colorbar_ticks=[0:20:100;],
     colorbar_title=L"t/T",
     colorbar_titlefontsize=22,
    #  colorbar_tickfontsize=2,
     labelfontsize=22,
     tickfontsize=12,
     margin=0.3 * Plots.cm)
title!(L"L=%$(L), n=%$(n), T=%$(int2label(t_meas))", titlefontsize=20)
xlabel!(L"\mathrm{site} \enspace i")
ylabel!(L"\left\langle n_{\mathtt{b},i} \right\rangle")
savefig("pics/bs/bs_expval_longword_L$(L)_n$(n)_w$(w)_$(int2str(t_meas)).pdf")


# <n_b>^2 over time
L = 120
n = L ÷ 10
w = 1
# data = deserialize("data/bs/bs_one_runs/bs_L$(L)_n$(n)_one.dat")
# data = deserialize("data/bs/bs_long_runs/bs_L$(L)_long_10^4.dat")
data = deserialize("data/bs/bs_b_wave_a_random/bs_L$(L)_n$(n)_w$(w)_10^3.dat")
# data = deserialize("data/bs/bs_random_words/thermalization/bs_L$(L)_10^3.dat")
# data = deserialize("data/bs/bs_n$(n)_w$(w)/bs_L$(L)_n$(n)_w$(w)_10^3.dat")
e = 10
exp_vals = data["nb2_all"][e]
L = data["L"]
n = data["n"]
# w = data["n_waves"]
t_therm = data["t_therm"]
t_meas = data["t_meas"]
tt = data["tt"][e] / t_meas
y_data = sum.((x -> x .^ 2).(exp_vals)) / L
roll_window = 30
smooth_y_data = rollmean(y_data, roll_window)
plot([1:length(y_data);], y_data, 
     legend=false, 
     lw=3,
     labelfontsize=22,
     tickfontsize=12,
     margin=0.3 * Plots.cm)
# tt = findfirst(x -> x ≤ y_data[1]*0.1, y_data)
# plot!([1:tt;], y_data[1]*(log2.((tt .- [1:tt;]) .+ 1) / log2(tt+1)), 
#       legend=false, 
#       lw=3)

# nt = 1/2*log2.(2^(2n) * (1 .- [1:tt;]/tt) .+ [1:tt;]/tt)
# nb2t = 4*nt.^3 / (L*(L/4-1))
# plot!([1:tt;], nb2t * y_data[1]/nb2t[1], 
#       legend=false, 
#       lw=3)
# plot!([tt]; seriestype = :vline)
# plot!([0.0025]; color=:orange, seriestype = :hline)
# plot!([y_data[1]*0.0175]; color=:goldenrod, seriestype = :hline)
# plot!([y_data[1]*0.10]; color=:lightgreen, seriestype = :hline)
# plot!([0.08]; color=:goldenrod, seriestype = :hline, lw=3)
plot!([1:length(smooth_y_data);] .+ (roll_window÷2), smooth_y_data, legend=false, lw=2, color="red")
xlabel!(L"t/T")
ylabel!(L"\overline{\left\langle n_{\mathtt{b}} \right\rangle^2}")
title!(L"L=%$(L), n=%$(n), T=%$(int2label(t_meas))", titlefontsize=20)
# title!(L"L=%$(L), T=%$(int2label(t_meas))," * "   initial state: ee...e", titlefontsize=16)
# title!(L"L=%$(L), t_\mathrm{cycle}=%$(int2label(t_meas)), \mathrm{initial \enspace state: random \enspace} w \sim e", titlefontsize=14)
# savefig("pics/bs/bs_nb2_eee_L$(L)_long_$(int2str(t_meas))")
savefig("pics/bs/bs_nb2_L$(L)_n$(n)_$(int2str(t_meas)).pdf")


# Thermalization time (1/10) vs system size L (for fixed n)
means = Float64[]
errors = Float64[]
L_list = [24:26; 28:2:40; 45:5:120;]
n = 5
for L in L_list
    data = deserialize("data/bs/bs_n$(n)_w1/bs_L$(L)_n$(n)_w1_30.dat")
    # println(length(data["exp_vals_all"][1]))
    tt = data["tt"]
    # From the long run calculation
    # if L == 50
    #     tt = [474600000]
    # end
    tt = tt[tt.>0]
    println("L: $(L), experiments thermalized: $(length(tt))")
    push!(means, mean(tt))
    push!(errors, std(tt))
end
# linear_fit = Polynomials.fit(L_list[3:end], log.(means[3:end]), 1)
# println(linear_fit.coeffs)
plot(L_list, means/10^4, yerror=errors/10^4,
    #  label=L"%$(round(exp(linear_fit.coeffs[1]); digits=2)) \cdot \exp(%$(round(linear_fit.coeffs[2]; digits=2))L)", 
    # yaxis=:log,
    # yticks=[10^4, 10^5, 10^6, 10^7, 10^8], 
    xticks=[20:10:120;],
    seriestype=:scatter,
    markersize=6,
    markerstrokewidth=2,
    color=:orange,
    label=false,
    labelfontsize=22,
    tickfontsize=16,
    margin=0.3 * Plots.cm)
# plot!(L_list[3:end], exp.(linear_fit.(L_list[3:end]))/10^6, 
#     # label=L"0.67 \cdot 2^{2n}",
#     color=:blue,
#     legend=:topright,
#     legend_font=22,
#     lw=2, 
#     labelfontsize=22,
#     tickfontsize=16)
xlabel!(L"$L$")
ylabel!(L"t_\mathrm{th} / 10^4")
title!(L"n=%$(n)", titlefontsize=22)
savefig("pics/bs/bs_therm_time_n$(n).pdf")


# Thermalization time (1/10) vs system size L (n = L/10)
means = Float64[]
errors = Float64[]
for L in 60:10:130
    data = deserialize("data/bs/bs_L_div_n_10/bs_L$(L)_n$(L÷10).dat")
    tt = data["tt"]
    push!(means, mean(tt))
    push!(errors, std(tt))
end
linear_fit = Polynomials.fit([60:10:130;], log.(means), 1)
println(linear_fit.coeffs)
plot([60:10:130;], (x -> 0.67 * 2^(2x/10)).([60:10:130;]), 
      label=L"0.67 \cdot 2^{2n}",
      color=:blue,
      legend=:bottomright,
      legend_font=22,
      lw=2, 
      labelfontsize=22,
      tickfontsize=16)
plot!([60:10:130;], means, yerror=errors,
    seriestype=:scatter,
    markersize=6,
    markerstrokewidth=2,
    color=:red,
    yaxis=:log, 
    yticks=[10^3, 10^4, 10^5, 10^6, 10^7, 10^8], 
    xticks=[60:10:130;],
    label=false)
xlabel!(L"$L$")
ylabel!(L"$t_\mathrm{th}$")
title!(L"n=L/10", titlefontsize=22)
savefig("pics/bs/bs_therm_time.pdf")


# Thermalization time vs number of waves w (for fixed L and n)
means = Float64[]
errors = Float64[]
L = 180
n = 8
w_list = [2, 3, 4]
therm_fraction = 0.10
roll_window = 30
for w in w_list
    data = deserialize("data/bs/bs_L$(L)_n$(n)/bs_L$(L)_n$(n)_w$(w)_300.dat")
    t_meas = data["t_meas"]
    nb2 = (exp_vals -> sum.((x -> x .^ 2).(exp_vals)) / L).(data["exp_vals_all"])
    smooth_nb2 = rollmean.(nb2, roll_window)
    
    tt_smooth = (y -> findfirst(x -> x ≤ y[1]*therm_fraction, y)).(smooth_nb2)
    tt_smooth = tt_smooth[.!(isnothing.(tt_smooth))]*t_meas
    println("w: $(w), experiments thermalized: $(length(tt_smooth))")
    push!(means, mean(tt_smooth))
    push!(errors, std(tt_smooth))
end
# linear_fit = Polynomials.fit(L_list, log.(means), 1)
# println(linear_fit.coeffs)
plot(w_list, means/10^4, yerror=errors/10^4,
    seriestype=:scatter,
    markersize=6,
    markerstrokewidth=2,
    color=:violet,
    # yticks=[10^4, 10^5, 10^6, 10^7, 10^8], 
    xticks=w_list,
    legend=false,
    legend_font=14,
    lw=2,
    labelfontsize=22,
    tickfontsize=16
    )
xlabel!(L"$w$")
ylabel!(L"$t_\mathrm{th} / 10^4$")
title!(L"L=%$(L), n=%$(n)", titlefontsize=22)
savefig("pics/bs/bs_therm_time_vs_nwaves_L$(L)_n$(n).pdf")


# Thermalization time vs n and number of waves w (for fixed L and fixed n*w, i.e., total amount of b's)
means = Float64[]
errors = Float64[]
L = 1700
n_list = [6,8,9,10]
w_list = 360 .÷ n_list
t_meas_list = [30, 10^3, 3*10^3, 3*10^4]
therm_fraction = 0.10
roll_window = 30
for (n, w, t_meas) in zip(n_list, w_list, t_meas_list)
    data = deserialize("data/bs/bs_contrast/bs_L$(L)_n$(n)_w$(w)_$(int2str(t_meas)).dat")
    nb2 = (exp_vals -> sum.((x -> x .^ 2).(exp_vals)) / L).(data["exp_vals_all"])
    smooth_nb2 = rollmean.(nb2, roll_window)
    
    tt_smooth = (y -> findfirst(x -> x ≤ y[1]*therm_fraction, y)).(smooth_nb2)
    tt_smooth = tt_smooth[.!(isnothing.(tt_smooth))]*t_meas
    println("L: $(L), n: $(n), w: $(w), t_meas: $(t_meas), experiments thermalized: $(length(tt_smooth))")
    push!(means, mean(tt_smooth))
    push!(errors, std(tt_smooth))
end
linear_fit = Polynomials.fit(n_list, log.(means), 1)
println(linear_fit.coeffs)
plot(n_list, 0.01 * 2 .^ (3.257*n_list),
      lw=2,
      c=:blue,
      label=L"0.01 \cdot 2^{3.257 n}"
      )
plot!(n_list, means, yerror=errors,
    yaxis=:log,
    seriestype=:scatter,
    c=:red,
    markersize=6,
    markerstrokewidth=2,
    yticks=[10^4, 10^5, 10^6, 10^7, 10^8], 
    xticks=[6:12;],
    label=false,
    legendfontsize=16,
    legend=:topleft,
    lw=2,
    labelfontsize=22,
    tickfontsize=16)
xlabel!(L"$n$")
ylabel!(L"$t_\mathrm{th}$")
title!(L"L=%$(L), w\cdot n = 360", titlefontsize=22)
savefig("pics/bs/bs_therm_time_fixed_contrast_L1700.pdf")



# Expansion length vs n fow w_large(n)
n_list = [5:10;]
EL = [24, 29, 34, 40, 44, 50]
linear_fit = Polynomials.fit(n_list, EL, 1)
println(linear_fit.coeffs)
plot(n_list, (x -> 5.2x-2).(n_list), 
      label=L"5.2n - 2",
      color=:blue,
      legend=:bottomright,
      legend_font=26,
      lw=4, 
      labelfontsize=26,
      tickfontsize=20)
plot!(n_list, EL,
    seriestype=:scatter,
    markersize=8,
    markerstrokewidth=4,
    color=:red,
    # yaxis=:log, 
    # yticks=[10^3, 10^4, 10^5, 10^6, 10^7, 10^8], 
    xticks=n_list,
    label=false)
xlabel!(L"$n$")
ylabel!(L"$\mathrm{EL}$")
# title!(L"n=L/10", titlefontsize=22)
savefig("pics/bs/bs_expansion_length.pdf")



# Animation
L = 130
data = deserialize("data/bs/bs_one_runs/bs_L$(L)_n$(L÷10)_one.dat")
exp_vals = data["exp_vals_all"][1]
nb2 = sum.((x -> x .^ 2).(exp_vals)) / L
n = data["n"]
t_therm = data["t_therm"]
t_meas = data["t_meas"]
anim = @animate for t in 1:100
    p1 = plot([1:L;], exp_vals[t], legend=false, lw=4, labelfontsize=32, xlabel=L"\mathrm{site} \enspace i", ylabel=L"\left\langle n_{\mathtt{b},i} \right\rangle", ylims=(-0.5, 0.55),
        title=L"t/T=%$(t),  T=%$(int2label(t_meas))", titlefont=24, tickfontsize=20)
    p2 = plot([1:t;], nb2[1:t], legend=false, lw=4, labelfontsize=32, xlabel=L"t/T", ylabel=L"\overline{\left\langle n_{\mathtt{b}} \right\rangle^2}",
        ylims=(0.0, 0.11), xlims=(0, 100), title=L"t/T=%$(t),  T=%$(int2label(t_meas))", titlefont=24, tickfontsize=20)
    plot(p1, p2, windowsize=(2000, 800), margin=1.5 * Plots.cm)
end
gif(anim, "pics/bs/anim_$(L).gif", fps=7)


# <nb>^2 over normalized time for different L, for a single experiment
Plots.CURRENT_PLOT.nullableplot = nothing
color_palette = palette(:rainbow, 8)
i = 1
for L in 60:10:130
    data = deserialize("data/bs/bs_L$(L)_n$(L÷10)_one.dat")
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
savefig("pics/bs/bs_nb2_scaled_time_one")


# <nb>^2 over normalized time for different L, where <nb>^2 is averaged over many experiments
Plots.CURRENT_PLOT.nullableplot = nothing
color_palette = palette(:rainbow, 8)
i = 1
for L in [120]
    data = deserialize("data/bs/bs_L_div_n_10/bs_L$(L)_n$(L÷10).dat")
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
ylabel!(L"\overline{\left\langle n_{\mathtt{b}} \right\rangle^2}")
title!("averaged over runs")
xlims!((0.0, 1.05))
savefig("pics/bs/bs_nb2_scaled_time")


# <nb>^2 over normalized time for L=120, together with the analytic formula
Plots.CURRENT_PLOT.nullableplot = nothing
L = 120
data = deserialize("data/bs/bs_L_div_n_10/bs_L$(L)_n$(L÷10).dat")
data_one = deserialize("data/bs/bs_one_runs/bs_L$(L)_n$(L÷10)_one.dat")
t_meas = data["t_meas"]
experiments = data["experiments"]
exp_vals_all = data["exp_vals_all"]
nb2_all = (exp_vals -> sum.((x -> x .^ 2).(exp_vals)) / L).(exp_vals_all)
nb2_one = sum.((x -> x .^ 2).(data_one["exp_vals_all"][1])) / L
tc_all = (nb2 -> t_meas * findfirst(x -> x ≤ nb2[1] / 10, nb2)).(nb2_all)
tc_one = data_one["t_meas"] * findfirst(x -> x ≤ nb2_one[1] / 10, nb2_one)
t_all = [[1:length(nb2);] * t_meas / tc for (nb2, tc) in zip(nb2_all, tc_all)]
t_one = [1:length(nb2_one);] * data_one["t_meas"] / tc_one
lin_int_all = [linear_interpolation(t, nb) for (t, nb) in zip(t_all, nb2_all)]
lin_int_one = linear_interpolation(t_one, nb2_one)
t = [t_meas/minimum(tc_all):0.001:1;]
nb2_all_interpolated = [lin_int(t) for lin_int in lin_int_all]
nb2_all_mean = mean(nb2_all_interpolated)
nb2_all_errors = std(nb2_all_interpolated)
plot(t, nb2_all_mean .- nb2_all_errors, fillrange = nb2_all_mean .+ nb2_all_errors,
      fillalpha = 0.35,
      fillcolor=:red,
      color=false,
      label="standard deviation")
plot!(t_one, nb2_one, color=:grey, alpha=0.5, lw=2, label="one run")
plot!(t, nb2_all_mean,
      color=:blue, 
      lw=2, 
      legend=:topright,
      labelfontsize=24,
      tickfontsize=14,
      margin=0.3 * Plots.cm,
      legendfontsize=12,
      label="averaged over runs")

tt = 1
n = L ÷ 10
nt = @. n + 1/2*log2(1 - t/tt + 2^(-2*Float64(n))*t/tt)
nb2t = nt.^2
plot!(t, nb2t * nb2_all_mean[1]/nb2t[1], 
    label="analytic formula",
    color=:green,
    lw=2,
    linestyle=:dash)
xlabel!(L"t/t_\mathrm{th}")
ylabel!(L"\overline{\left\langle n_{\mathtt{b}} \right\rangle^2}")
# title!("averaged over runs")
xlims!((0.0, t_one[end]))
xticks!([0.0:0.25:1.75;])
savefig("pics/bs/bs_nb2_scaled_time.pdf")


# Plot the number of group elements vs L
n_group_elems = deserialize("data/bs/n_group_elems.dat")
linear_fit = Polynomials.fit([5:length(n_group_elems);], log.(3, n_group_elems[5:end]), 1)
coeffs = Float64.(linear_fit.coeffs)
println(linear_fit)
plot([1:length(n_group_elems);], n_group_elems, seriestype=:scatter, yaxis=:log,
      legend=:topleft, label="exact", yticks=[1, 10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8], markersize=6, markerstrokewidth=2)
plot!([1:length(n_group_elems);], (3).^(coeffs[1] .+ [1:length(n_group_elems);] / 2), 
      color=:red, lw=3, legend_font=18, labelfontsize=20, tickfontsize=16,
      label=L"%$(round(3^(coeffs[1]); digits=2)) \cdot 3^{L/2}")
xlabel!(L"$L$")
ylabel!(L"N_K")
savefig("pics/bs/bs_group_volume.pdf")


# Plot the number of words in the identity sector vs L
n_id_words = deserialize("data/bs/n_id_words.dat")
Ls = [1:length(n_id_words);]
n_id_fraction = n_id_words ./ (5 .^ Ls)
# linear_fit = Polynomials.fit(Ls, log.(n_id_words), 1)
# coeffs = Float64.(linear_fit.coeffs)
# println(coeffs)
# plot(Ls, n_id_words, seriestype=:scatter, legend=:topleft, yaxis=:log,
#       label="exact", yticks=[1, 10^2, 10^4, 10^6, 10^8, 10^10, 10^12, 10^14, 10^16])
# plot!(Ls, exp.(linear_fit.(Ls)), 
#       color=:red, lw=2, legend_font=12, labelfontsize=16, tickfontsize=10,
#       label=L"%$(round(exp(coeffs[1]); digits=2)) \cdot \exp(%$(round(coeffs[2]; digits=2))L)",
#       margin=0.3 * Plots.cm)
linear_fit = Polynomials.fit(Ls .^ (1/3), log.(n_id_fraction), 1)
coeffs = Float64.(linear_fit.coeffs)
println(coeffs)
plot(Ls, n_id_fraction, seriestype=:scatter, legend=:topright, yaxis=:log,
      label="exact", yticks=[1, 10^(-1), 10^(-2), 10^(-2), 10^(-3)])
plot!(Ls, exp.(linear_fit.(Ls .^ (1/3))), 
      color=:red, lw=2, legend_font=12, labelfontsize=16, tickfontsize=10,
      label=L"%$(round(exp(coeffs[1]); digits=2)) \cdot \exp(%$(round(coeffs[2]; digits=2))L^{1/3})",
      margin=0.3 * Plots.cm)
xlabel!(L"$L$")
ylabel!("fraction of words \n in the identity sector")
savefig("pics/bs/n_id_fraction")


# Distribution of <nb>^2 averaged over 10^3 time steps starting from random words in the identity sector
L = 100
N = 10^6
t_meas = 10^3
Plots.CURRENT_PLOT.nullableplot = nothing
# colors = palette(:thermal, length(L_list))
nb2 = deserialize("data/bs/bs_random_words/identity_sector/bs_nb2_random_words_L$(L)_N$(int2str(N))_$(int2str(t_meas)).dat")
h = StatsBase.fit(Histogram, nb2, 0:0.00025:0.03)
edg = collect(h.edges[1])
x_data = (edg[1:end-1] .+ edg[2:end])/2
y_data = h.weights
plot(nb2, bins=edg, label=false, legend=:topright, minorgrid=true, legend_font=10, seriestype=:stephist, color=:black, 
        fill=true, fillcolor=:gray, fillalpha=0.75, labelfontsize=16, tickfontsize=12, margin=0.3 * Plots.cm)

# plot!([0.075]; seriestype = :vline)
# xlims!((0.0, 0.03))
xlabel!(L"\overline{\left\langle n_{\mathtt{b},i} \right\rangle^2}")
ylabel!("number of random words")
title!(L"\mathrm{identity} \enspace \mathrm{sector}, L=%$(L), T=%$(int2str(t_meas))", titlefontsize=18)
savefig("pics/bs/bs_nb2_distribution_L$(L)_$(int2str(t_meas))")


# Distribution of <nb>^2 averaged over 10^3 time steps starting from random words in the 0 b-charge sector
L = 100
N = 10^5
t_meas = 3*10^3
Plots.CURRENT_PLOT.nullableplot = nothing
# colors = palette(:thermal, length(L_list))
nb2 = deserialize("data/bs/bs_random_words/0b_sector/bs_nb2_random_words_L$(L)_N$(int2str(N))_$(int2str(t_meas)).dat")
h = StatsBase.fit(Histogram, nb2, 0:0.00025:0.03)
edg = collect(h.edges[1])
x_data = (edg[1:end-1] .+ edg[2:end])/2
y_data = h.weights
plot(nb2, bins=edg, label=false, legend=:topright, minorgrid=true, legend_font=10, seriestype=:stephist, color=:black, 
        fill=true, fillcolor=:gray, fillalpha=0.75, labelfontsize=16, tickfontsize=12, margin=0.3 * Plots.cm)
# plot!([0.075]; seriestype = :vline)
# xlims!((0.0, 0.03))
xlabel!(L"\overline{\left\langle n_{\mathtt{b},i} \right\rangle^2}")
ylabel!("number of random words")
# title!(L"\mathrm{identity} \enspace \mathrm{sector}, L=%$(L), T=%$(int2str(t_meas))", titlefontsize=18)
title!(L"0 \enspace \mathrm{b}\mathrm{-}\mathrm{charge} \enspace \mathrm{sector}, L=%$(L), T=%$(int2label(t_meas))", titlefontsize=18)
savefig("pics/bs/bs_nb2_distribution_0bcharge_L$(L)_$(int2str(t_meas)).pdf")


# nb averaged over random words in the identity sector
L = 100
nb = deserialize("data\\bs\\bs_random_words\\nb_id_sector\\bs_nb_random_words_L$(L)_N10^5.dat")
plot([1:L;], nb, legend=false, minorgrid=true, labelfontsize=18, tickfontsize=16, margin=0.3 * Plots.cm, lw=3)
# title!(L"L=%$(L), \mathrm{averaged} \enspace \mathrm{over} \enspace 10^5 \enspace \mathrm{states}", titlefontsize=18)
ylabel!(L"\left\langle n_{\mathtt{b},i} \right\rangle")
xlabel!(L"site $i$")
savefig("pics/bs/bs_nb_L$(L).pdf")



# Thermalization time vs L (n=L/10) for waves of b with random stuff in between
medians = Float64[]
means = Float64[]
L_list = [50:10:130;]
t_meas_list = [1,1,3,10,30,100,300,10^3,3*10^3]
therm_fraction = 0.75
roll_window = 100
for (L, t_meas) in zip(L_list, t_meas_list)
    data = deserialize("data/bs/bs_b_wave_a_random/bs_L$(L)_n$(L÷10)_w1_$(int2str(t_meas)).dat")
    nb2_all = data["nb2_all"]
    smooth_nb2 = rollmean.(nb2_all, roll_window)
    
    tt_smooth = (y -> findfirst(x -> x ≤ y[1]*therm_fraction, y)).(smooth_nb2)
    tt_smooth = tt_smooth[.!(isnothing.(tt_smooth))]*t_meas
    println("L: $(L), experiments thermalized: $(length(tt_smooth))")
    if length(tt_smooth) > 0
        push!(medians, median(tt_smooth))
        push!(means, mean(tt_smooth))
    else
        push!(medians, NaN)
        push!(means, NaN)
    end
end
# linear_fit = Polynomials.fit(L_list, log2.(medians), 1)
# println(linear_fit.coeffs)
plot(L_list,
    (x -> 9.314 * 2^(0.137x)).(L_list), 
    label=L"0.67 \cdot 2^{1.37n}",
    color=:blue,
    legend=:bottomright,
    legend_font=22,
    lw=2, 
    labelfontsize=22,
    tickfontsize=16)
plot!(L_list, medians, 
    # yerror=errors,
    seriestype=:scatter,
    markersize=6,
    markerstrokewidth=2,
    color=:lime,
    yaxis=:log,
    label=false,
    yticks=[10^4, 10^5, 10^6, 10^7, 10^8], 
    xticks=L_list,
    labelfontsize=22,
    tickfontsize=16
    )
xlabel!(L"$L$")
ylabel!(L"$t_\mathrm{th}$")
# title!(L"L=%$(L), n=L/10, ", titlefontsize=22)
savefig("pics/bs/bs_t_th_b_wave_a_random.pdf")



# nb2 (n=L/10) for waves of b with random stuff in between
L = 100
t_meas = 100
therm_fraction = 0.75
roll_window = 100
data = deserialize("data/bs/bs_b_wave_a_random/bs_L$(L)_n$(L÷10)_w1_$(int2str(t_meas))_new.dat")
# experiments = [30,24,35,11,2,3,4,6,14,15,16,26]
experiments = [2,36,18,50,12,4,23,44,1,3,14,21,30]
exp_flat = [1,3,14,21,30]
exp_drop = [2,4,12,18,23, 36, 44, 50]
nb2_all = data["nb2_all"][experiments]
smooth_nb2_all = rollmean.(nb2_all, roll_window)
println(length.(nb2_all))
Plots.CURRENT_PLOT.nullableplot = nothing
red_palette = palette(:heat,256)
blue_palette = palette([:turquoise, :blue],256)
i = 1
for (nb2, smooth_nb2, experiment) in zip(nb2_all, smooth_nb2_all, experiments)
    if experiment in exp_drop
        color = red_palette[32 + 28i]
        i += 1
    elseif experiment in exp_flat
        color = :blue
    end
    plot!([1:length(smooth_nb2)], smooth_nb2, 
        legend=false, 
        lw=3,
        labelfontsize=22,
        tickfontsize=12,
        margin=0.3 * Plots.cm,
        color=color
        )
end
xlabel!(L"t/T")
ylabel!(L"\overline{\left\langle n_{\mathtt{b}} \right\rangle^2}")
# title!(L"L=%$(L), n=%$(n), T=%$(int2label(t_meas))", titlefontsize=20)
savefig("pics/bs/bs_nb2_b_wave_a_random.pdf")


