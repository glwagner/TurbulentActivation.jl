# Overlay the online chamber's activation curves with the reference replay through the SAM
# fields (same droplet physics, same counterfactual construction), with the spread across
# windows as bands.
#
# Usage: julia --project analysis/compare_curves.jl online_activation.jld2 reference_replay.jld2 [output.png]
using JLD2, CairoMakie, Statistics

online = load(ARGS[1])["results"]
reference = load(ARGS[2])["results"]
output = length(ARGS) ≥ 3 ? ARGS[3] : "comparison.png"

curves(results) = begin
    targets = 100 .* results["targets"]
    windows = [k for k in keys(results) if startswith(k, "window_")]
    stack(kind) = reduce(hcat, [results[w][kind] for w in windows])
    (; targets, windows, Dict(kind => (mean(stack(kind), dims=2)[:], length(windows) > 1 ? std(stack(kind), dims=2)[:] : zeros(length(targets))))
                              for kind in ("fluctuating", "uniform", "instantaneous"))...)
end
o = curves(online); r = curves(reference)
@info "comparison" online_windows=length(o.windows) reference_windows=length(r.windows)

fig = Figure(size=(760, 500), fontsize=16)
ax = Axis(fig[1, 1]; xlabel="mean supersaturation (%)", ylabel="activated fraction after 60 s",
          title="Online Breeze chamber (solid) vs. replay through the SAM fields (dashed)")
colors = Dict("fluctuating" => :dodgerblue3, "uniform" => :gray40, "instantaneous" => :darkorange2)
for kind in ("uniform", "instantaneous", "fluctuating")
    m, s = getproperty(o, Symbol(kind)); band!(ax, o.targets, m .- s, m .+ s; color=(colors[kind], 0.15))
    lines!(ax, o.targets, m; color=colors[kind], linewidth=2.5, label="$kind, Breeze online")
    m, s = getproperty(r, Symbol(kind)); band!(ax, r.targets, m .- s, m .+ s; color=(colors[kind], 0.15))
    lines!(ax, r.targets, m; color=colors[kind], linewidth=2.5, linestyle=:dash, label="$kind, SAM replay")
end
vlines!(ax, [0.0]; color=:black, linestyle=:dot)
axislegend(ax; position=:lt, nbanks=1, labelsize=12)
save(output, fig)
@info "saved" output
