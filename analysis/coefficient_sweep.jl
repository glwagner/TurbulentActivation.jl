# The wall transfer coefficient sweep: fluctuating and instantaneous activation curves of the online
# chamber for several bulk coefficients against the replay through the reference SAM fields.
#
# Usage: julia --project analysis/coefficient_sweep.jl reference_replay.jld2 output.png C₁=file₁.jld2 C₂=file₂.jld2 ...
using JLD2, CairoMakie, Statistics

reference = load(ARGS[1])["results"]
output = ARGS[2]
runs = [(label = rsplit(a, "="; limit=2)[1], results = load(rsplit(a, "="; limit=2)[2])["results"]) for a in ARGS[3:end]]

function mean_curve(results, kind)
    windows = sort([k for k in keys(results) if startswith(k, "window_")])
    A = reduce(hcat, [results[w][kind] for w in windows])
    return 100 .* results["targets"], mean(A, dims=2)[:]
end

fig = Figure(size=(1100, 460), fontsize=15)
colors = Makie.wong_colors()
for (i, kind) in enumerate(("fluctuating", "instantaneous"))
    ax = Axis(fig[1, i]; xlabel="mean supersaturation (%)", ylabel="activated fraction after 60 s",
              title=kind * " replicas", limits=((-4, 2), (-0.02, 1.02)))
    x, y = mean_curve(reference, kind)
    lines!(ax, x, y; color=:black, linestyle=:dash, linewidth=3, label="replay through the SAM fields")
    for (j, run) in enumerate(runs)
        x, y = mean_curve(run.results, kind)
        lines!(ax, x, y; color=colors[j], linewidth=2.5, label="Breeze chamber, C = " * run.label)
    end
    vlines!(ax, [0.0]; color=:gray50, linestyle=:dot)
    i == 1 && axislegend(ax; position=:lt, labelsize=12)
end
save(output, fig)
@info "saved" output
