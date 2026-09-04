# Eulerian statistics of the reference SAM fields: domain mean and standard deviation of the
# supersaturation (Anderson's definition, Magnus e_sat and e = p q / (0.622 + q)), temperature,
# vapor, and condensate, per snapshot and averaged over the sampled snapshots, with the
# near-wall cells excluded from an "interior" variant.
#
# Usage: julia --project analysis/reference_statistics.jl [stride] [output.jld2]
using JLD2, Statistics, Printf
include(joinpath(@__DIR__, "..", "reference", "read_bin3d.jl"))
using .SAMBin3D

dir = joinpath(homedir(), "anderson_reference", "LES", "LES")
files = sort(filter(f -> endswith(f, ".bin3D"), readdir(dir; join=true)))
stride = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 20
output = length(ARGS) ≥ 2 ? ARGS[2] : joinpath(@__DIR__, "..", "runs", "reference_statistics.jld2")
selected = files[1:stride:end]
@info "Reference statistics" nfiles=length(files) sampled=length(selected) stride

esat(T) = 610.94 * exp(17.625 * (T - 273.15) / (T - 30.11))          # Pa, T in K
function supersaturation(T, qv, p)
    e = p * qv / (0.622 + qv)
    return e / esat(T) - 1
end

rows = []
for (n, file) in enumerate(selected)
    s = read_bin3d(file)
    T = Float64.(s.TABS); qv = Float64.(s.QV) ./ 1000; qn = Float64.(s.QN) ./ 1000
    p = Float64.(s.p) .* 100
    S = similar(T)
    for k in axes(T, 3), j in axes(T, 2), i in axes(T, 1)
        S[i, j, k] = supersaturation(T[i, j, k], qv[i, j, k], p[k])
    end
    interior = @view S[3:end-2, 3:end-2, 3:end-2]
    row = (; step = parse(Int, match(r"_(\d+)\.bin3D", file).captures[1]), time = s.time,
             S_mean = mean(S), S_std = std(S), S_interior_mean = mean(interior), S_interior_std = std(interior),
             S_profile = [mean(S[:, :, k]) for k in axes(S, 3)],
             T_mean = mean(T), T_std = std(T), qv_mean = mean(qv), qv_std = std(qv),
             qn_mean = mean(qn), qn_max = maximum(qn), cloudy_fraction = mean(qn .> 1e-6))
    push!(rows, row)
    n % 10 == 1 && @printf("%s: ⟨S⟩ = %+.3f %%  σ(S) = %.3f %%  interior σ = %.3f %%  ⟨T⟩ = %.2f K  σ(T) = %.3f K  ⟨qn⟩ = %.3f g/kg\n",
                          basename(file), 100row.S_mean, 100row.S_std, 100row.S_interior_std, row.T_mean, row.T_std, 1000row.qn_mean)
end

summary = (; S_mean = mean(r.S_mean for r in rows), S_std = mean(r.S_std for r in rows),
             S_interior_mean = mean(r.S_interior_mean for r in rows), S_interior_std = mean(r.S_interior_std for r in rows),
             T_std = mean(r.T_std for r in rows), qv_std = mean(r.qv_std for r in rows), qn_mean = mean(r.qn_mean for r in rows),
             S_profile = mean(hcat((r.S_profile for r in rows)...), dims=2)[:], z = read_bin3d(first(selected)).z)
@printf("\nreference over %d snapshots: ⟨S⟩ = %+.3f %%  σ(S) = %.3f %%  (interior: %+.3f %%, %.3f %%)  σ(T) = %.3f K  σ(qv) = %.3f g/kg  ⟨qn⟩ = %.3f g/kg\n",
        length(rows), 100summary.S_mean, 100summary.S_std, 100summary.S_interior_mean, 100summary.S_interior_std, summary.T_std, 1000summary.qv_std, 1000summary.qn_mean)
jldsave(output; rows, summary)
@info "saved" output
