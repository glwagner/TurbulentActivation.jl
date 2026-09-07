# Compare Breeze's κ-Köhler droplet growth with Anderson's own code on identical supersaturation
# traces: `analysis/export_traces.jl` writes the traces, `analysis/anderson_growth.py` integrates
# them with `lfierce2/LagrangianDroplets`, and this script integrates the same traces with
# `Breeze.LagrangianMicrophysics.implicit_growth_step` and compares diameters, activation
# fractions, and the critical point. This is the link between our replay and Anderson's paper:
# the replay is his experiment run with our droplet physics, and this measures that substitution.
#
# Usage: julia --project analysis/growth_comparison.jl [tracedir] [figure.png]
using Breeze, CairoMakie, DelimitedFiles, Printf, Statistics
using Breeze.LagrangianMicrophysics: implicit_growth_step, equilibrium_diameter, critical_diameter, critical_supersaturation

tracedir = length(ARGS) ≥ 1 ? ARGS[1] : "runs/growth_comparison"
figure = length(ARGS) ≥ 2 ? ARGS[2] : "figures/growth_comparison.png"

FT = Float64
constants = ThermodynamicConstants(FT)
Dᵈ, κ, T, p = FT(130e-9), FT(1), FT(287.6), FT(1e5)
parameters = DropletDynamics(FT; substeps=25)
Dᶜ = critical_diameter(Dᵈ, κ, T, constants)
𝒮ᶜ = critical_supersaturation(Dᵈ, κ, T, constants)
@printf("Breeze: critical diameter %.4e m, critical supersaturation %.5f %%\n", Dᶜ, 100𝒮ᶜ)

grow(D², 𝒮, Δt) = implicit_growth_step(D², 𝒮, T, p, Dᵈ, κ, Δt, parameters, constants)

fig = Figure(size=(1000, 380), fontsize=15)
targets = String[]
rows = []
axes = Axis[]
for (i, path) in enumerate(sort(filter(f -> startswith(basename(f), "traces_"), readdir(tracedir, join=true))))
    label = replace(basename(path), "traces_" => "", ".csv" => "")
    A = readdlm(path, ',', FT)
    t, S = A[:, 1], A[:, 2:end]
    dpath = replace(path, "traces_" => "diameters_")
    isfile(dpath) || continue
    B = readdlm(dpath, ',', FT)
    Dₐ = B[:, 2:end]                       # Anderson's diameters
    n = size(S, 2)
    Dᵦ = zeros(FT, size(Dₐ))               # ours, from his initial size
    Dᵦᵒ = zeros(FT, size(Dₐ))              # ours, from our own initial size
    for m in 1:n
        D² = Dₐ[1, m]^2
        D²ᵒ = equilibrium_diameter(min(S[1, m], FT(0.9) * 𝒮ᶜ), Dᵈ, κ, T, constants)^2
        Dᵦ[1, m], Dᵦᵒ[1, m] = sqrt(D²), sqrt(D²ᵒ)
        for k in 2:length(t)
            Δt = t[k] - t[k-1]
            D² = grow(D², S[k, m], Δt)
            D²ᵒ = grow(D²ᵒ, S[k, m], Δt)
            Dᵦ[k, m], Dᵦᵒ[k, m] = sqrt(D²), sqrt(D²ᵒ)
        end
    end
    fₐ = count(≥(Dᶜ), Dₐ[end, :]) / n
    fᵦ = count(≥(Dᶜ), Dᵦ[end, :]) / n
    fᵦᵒ = count(≥(Dᶜ), Dᵦᵒ[end, :]) / n
    Δ = 100 * maximum(abs.(Dᵦ .- Dₐ) ./ max.(Dₐ, eps(FT)))
    push!(targets, label)
    push!(rows, (label, fₐ, fᵦ, fᵦᵒ, Δ, median(abs.(Dᵦ[end, :] .- Dₐ[end, :]) ./ Dₐ[end, :]) * 100))
    ax = Axis(fig[1, i]; xlabel="time in window (s)", ylabel=i == 1 ? "droplet diameter (μm)" : "",
              title="target mean 𝒮 = $label %", yscale=log10)
    for m in 1:min(n, 40)
        lines!(ax, t, 1e6 .* Dₐ[:, m]; color=(:black, 0.35), linewidth=1)
        lines!(ax, t, 1e6 .* Dᵦ[:, m]; color=(:dodgerblue3, 0.35), linewidth=1, linestyle=:dash)
    end
    hlines!(ax, [1e6 * Dᶜ]; color=:firebrick, linestyle=:dot)
    push!(axes, ax)
end
linkyaxes!(axes...)
for ax in axes[2:end]
    hideydecorations!(ax; grid=false)
end
Legend(fig[1, 4], [LineElement(color=:black), LineElement(color=:dodgerblue3, linestyle=:dash), LineElement(color=:firebrick, linestyle=:dot)],
       ["Anderson's code", "Breeze", "critical diameter"], framevisible=false)
save(figure, fig)

println("\n target   activated: Anderson   Breeze (his D₀)   Breeze (our D₀)   max|ΔD|/D   median|ΔD|/D at 60 s")
for (label, fₐ, fᵦ, fᵦᵒ, Δ, med) in rows
    @printf("  %6s        %.3f            %.3f             %.3f          %6.2f %%        %6.3f %%\n", label, fₐ, fᵦ, fᵦᵒ, Δ, med)
end
@info "saved" figure
