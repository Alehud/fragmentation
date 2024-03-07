include("../src/Fragmentation.jl")
using .Fragmentation
using IterTools
using SparseArrays
using Serialization
using Graphs
using Luxor
using Karnak
using Colors
using Plots
using LaTeXStrings

# H = 1/sqrt(2) * ComplexF64[1 1; 1 -1]
# T = ComplexF64[1 0; 0 exp(pi/4*im)]
# gate_set = Matrix{ComplexF64}[H, T]

X = ComplexF64[0 1; 1 0]
Y = ComplexF64[0 -im; im 0]
Z = ComplexF64[1 0; 0 -1]
Id = ComplexF64[1 0; 0 1]
θx = 0.01
θz = 0.01
Rx = exp(-im*X*θx)
Rz = exp(-im*Z*θz)
gate_set = Matrix{ComplexF64}[Rx, Id]

L_max = 10
alphas = Float64[]
betas = Float64[]
gammas = Float64[]
deltas = Float64[]
for L in 1:L_max
    for state in product(fill([0:1;], L)...)
        # println(state)
        U = ComplexF64[1 0; 0 1]
        for g in state
            U *= gate_set[g+1]
        end
        a11 = angle(U[1, 1])
        a12 = angle(U[1, 2])
        a22 = angle(U[2, 2])
        α = (a11 + a22)/2
        β = a12 - a11
        δ = a22 - a12
        γ = 2*acos(abs(U[1,1]))
        if γ < 0
            println(γ)
            flush(stdout)
        end
        push!(alphas, α)
        push!(betas, β)
        push!(gammas, γ)
        push!(deltas, δ)
    end
end

if gate_set[1] == H
    label = "HT"
elseif gate_set[1] == Rx
    label = "RxRz"
end
serialize("data/circuit/U_$(label)_params_L$(L_max).dat", Dict("α" => alphas, "β" => betas, "γ" => gammas, "δ" => deltas))

# d = deserialize("data/circuit/U_$(label)_params_L$(L_max).dat")
# alphas = d["α"]
# betas = d["β"]
# gammas = d["γ"]
# deltas = d["δ"]


plot(alphas, betas, seriestype=:scatter, color=:red, legend=false)
xlabel!(L"\alpha")
ylabel!(L"\beta")
savefig("pics/circuit/$(label)_alpha_beta_L$(L_max)")

plot(alphas, gammas, seriestype=:scatter, color=:red, legend=false)
xlabel!(L"\alpha")
ylabel!(L"\gamma")
savefig("pics/circuit/$(label)_alpha_gamma_L$(L_max)")

plot(alphas, deltas, seriestype=:scatter, color=:red, legend=false)
xlabel!(L"\alpha")
ylabel!(L"\delta")
savefig("pics/circuit/$(label)_alpha_delta_L$(L_max)")

plot(betas, gammas, seriestype=:scatter, color=:red, legend=false)
xlabel!(L"\beta")
ylabel!(L"\gamma")
savefig("pics/circuit/$(label)_beta_gamma_L$(L_max)")

plot(betas, deltas, seriestype=:scatter, color=:red, legend=false)
xlabel!(L"\beta")
ylabel!(L"\delta")
savefig("pics/circuit/$(label)_beta_delta_L$(L_max)")

plot(gammas, deltas, seriestype=:scatter, color=:red, legend=false)
xlabel!(L"\gamma")
ylabel!(L"\delta")
savefig("pics/circuit/$(label)_gamma_delta_L$(L_max)")


