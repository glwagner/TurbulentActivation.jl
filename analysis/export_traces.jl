# Export supersaturation traces from the reference replay for the cross-code growth comparison:
# the Lagrangian 𝒮(t) of the first `N` droplets of one window, shifted to Anderson's target
# means exactly as the fluctuating replicas are constructed, written as CSV for
# `analysis/anderson_growth.py` (his code) and read back by `analysis/growth_comparison.jl` (ours).
#
# Usage: julia --project analysis/export_traces.jl [N] [replay.jld2] [outdir]
using JLD2, Statistics, Printf, DelimitedFiles

N = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 200
replay = length(ARGS) ≥ 2 ? ARGS[2] : "runs/reference_replay.jld2"
outdir = length(ARGS) ≥ 3 ? ARGS[3] : "runs/growth_comparison"
mkpath(outdir)

f = load(replay)
series = f["series"]["window_1"]
S = series["S"]                      # (droplets, times)
t = series["sample_times"]
n = min(N, size(S, 1))
targets = (-0.02, -0.01, 0.0)
S̄ = mean(S)                          # the window mean removed in the fluctuating construction

@info "exporting traces" n length(t) mean_S=S̄ targets
for target in targets
    name = @sprintf("traces_%+05.1f.csv", 100target)
    A = hcat(t, permutedims(target .+ (S[1:n, :] .- S̄)))   # time in column 1, one droplet per column
    writedlm(joinpath(outdir, name), A, ',')
end
writedlm(joinpath(outdir, "times.csv"), t, ',')
println("EXPORT OK ", outdir)
