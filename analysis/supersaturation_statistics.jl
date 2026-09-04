# Lagrangian supersaturation statistics from the sampled droplets of the window experiment:
# the PDF of 𝒮, its Lagrangian autocorrelation and correlation time τₛ (Anderson et al. 2023
# report τₛ ≈ 7.5 s), and example trajectories of 𝒮 and D.
#
# Usage: julia --project analysis/supersaturation_statistics.jl <prefix>_activation.jld2 [output.png]
using JLD2, CairoMakie, Statistics

file = ARGS[1]
output = length(ARGS) ≥ 2 ? ARGS[2] : replace(file, "_activation.jld2" => "_supersaturation.png")
data = load(file)
series = data["series"]
windows = sort([k for k in keys(series) if haskey(series[k], "S")])
@info "Supersaturation statistics" file windows

# Autocorrelation of the fluctuation about the window-mean trajectory, averaged over droplets
function autocorrelation(S, Δt)                    # S is (droplets, times)
    S′ = S .- mean(S, dims=2)
    nt = size(S, 2)
    lags = 0:nt-1
    R = [mean(S′[:, 1:nt-l] .* S′[:, 1+l:nt]) for l in lags] ./ mean(S′ .^ 2)
    τ = Δt * (sum(R[R .> 0]) - 0.5)                 # integral time scale up to the first zero crossing
    return lags .* Δt, R, τ
end

fig = Figure(size=(1100, 700), fontsize=15)
ax1 = Axis(fig[1, 1]; xlabel="supersaturation (%)", ylabel="probability density", title="Lagrangian PDF of 𝒮")
ax2 = Axis(fig[1, 2]; xlabel="lag (s)", ylabel="autocorrelation", title="Lagrangian autocorrelation of 𝒮′")
ax3 = Axis(fig[2, 1:2]; xlabel="time in window (s)", ylabel="𝒮 (%)", title="Three droplet trajectories (first window)")
for (i, w) in enumerate(windows)
    S = series[w]["S"]'                               # (droplets, times)
    t = series[w]["sample_times"]
    Δt = length(t) > 1 ? t[2] - t[1] : 0.5
    hist!(ax1, 100 .* vec(S); bins=80, normalization=:pdf, color=(:dodgerblue3, 0.25))
    lags, R, τ = autocorrelation(S, Δt)
    lines!(ax2, lags, R; label="$(w): τₛ = $(round(τ, digits=1)) s, σ = $(round(100std(S), digits=2)) %")
    @info w mean_S=mean(S) std_S=std(S) τₛ=τ
    if i == 1
        for n in 1:3
            lines!(ax3, t, 100 .* S[n, :])
        end
        hlines!(ax3, [0.0]; color=:black, linestyle=:dot)
    end
end
axislegend(ax2; position=:rt)
save(output, fig)
@info "saved" output
