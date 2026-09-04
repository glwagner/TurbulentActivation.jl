# Activated fraction after 60 s versus target mean supersaturation, for the fluctuating,
# uniform, and instantaneous replicas, averaged over the sampling windows with the spread
# across windows as the uncertainty (Anderson et al. 2023, Fig. 3).
#
# Usage: julia --project analysis/activation_curve.jl <prefix>_activation.jld2 [output.png]
using JLD2, CairoMakie, Statistics

file = ARGS[1]
output = length(ARGS) ≥ 2 ? ARGS[2] : replace(file, ".jld2" => ".png")

data = load(file)
results = data["results"]
targets = 100 .* results["targets"]                       # percent
windows = [k for k in keys(results) if startswith(k, "window_")]
@info "Activation curve" file windows=length(windows)

stack(kind) = reduce(hcat, [results[w][kind] for w in windows])   # (targets, windows)
curves = Dict(kind => (mean(stack(kind), dims=2)[:], length(windows) > 1 ? std(stack(kind), dims=2)[:] : zeros(length(targets)))
              for kind in ("fluctuating", "uniform", "instantaneous"))

fig = Figure(size=(720, 480), fontsize=16)
ax = Axis(fig[1, 1]; xlabel="mean supersaturation (%)", ylabel="activated fraction after 60 s",
          title="Pi Chamber, $(length(windows)) window(s), N = $(results["N"]), grid $(results["size"])")
colors = Dict("fluctuating" => :dodgerblue3, "uniform" => :gray40, "instantaneous" => :darkorange2)
labels = Dict("fluctuating" => "fluctuating 𝒮", "uniform" => "uniform 𝒮", "instantaneous" => "instantaneous (τ = 0)")
for kind in ("uniform", "instantaneous", "fluctuating")
    m, s = curves[kind]
    band!(ax, targets, m .- s, m .+ s; color=(colors[kind], 0.2))
    scatterlines!(ax, targets, m; color=colors[kind], label=labels[kind], markersize=10)
end
vlines!(ax, [0.0]; color=:black, linestyle=:dot)
axislegend(ax; position=:lt)
save(output, fig)
@info "saved" output
