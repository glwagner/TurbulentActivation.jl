# Two summary figures: vertical profiles of the default chamber against the reference SAM fields,
# and a calibration map of every chamber run on the (Lagrangian spread, activation enhancement) plane.
# Usage: julia --project analysis/results_figures.jl
using JLD2, CairoMakie, Statistics, Printf
include(joinpath(@__DIR__, "..", "reference", "read_bin3d.jl"))
using .SAMBin3D

##### Profiles
d = jldopen("runs/weno5_loglaw1.25mm_profiles.jld2")
its = sort(parse.(Int, filter(k -> k != "serialized", keys(d["timeseries/T"]))))
last = its[end-36:end]                                     # the three 60 s windows (5 s cadence)
prof(k) = mean([vec(d["timeseries/$k/$i"]) for i in last])
T̄ = prof("T"); ℋ̄ = prof("ℋ"); ww = prof("ww")
keep = T̄ .> 0
T̄, ℋ̄, ww = T̄[keep], ℋ̄[keep], ww[keep]
close(d)
Nz = length(T̄); zc = ((1:Nz) .- 0.5) ./ Nz

dir = joinpath(homedir(), "anderson_reference", "LES", "LES")
Tr = zeros(32); Sr = zeros(32); wr = zeros(32); n = 0
esat(T) = 610.94 * exp(17.625 * (T - 273.15) / (T - 30.11))
for step in 90000:15000:180000
    s = read_bin3d(joinpath(dir, @sprintf("PiChamber_huji_19K_trj_32_%010d.bin3D", step)))
    q = s.QV ./ 1000; S = @. 1e5 * q / (0.622 + q) / esat(s.TABS) - 1
    for k in 1:32
        Tr[k] += mean(s.TABS[:, :, k]); Sr[k] += mean(S[:, :, k]); wr[k] += mean(s.W[:, :, k] .^ 2)
    end
    global n += 1
end
Tr ./= n; Sr ./= n; wr = sqrt.(wr ./ n); zr = ((1:32) .- 0.5) ./ 32

fig = Figure(size=(1000, 420), fontsize=15)
ax1 = Axis(fig[1, 1]; xlabel="mean supersaturation (%)", ylabel="z / H", title="⟨𝒮⟩(z)")
lines!(ax1, 100 .* (ℋ̄ .- 1), zc; linewidth=2.5, label="Breeze chamber")
lines!(ax1, 100 .* Sr, zr; linewidth=2.5, linestyle=:dash, color=:black, label="SAM reference")
vlines!(ax1, [0.0]; color=:gray60, linestyle=:dot)
axislegend(ax1; position=:rb, labelsize=12)
ax2 = Axis(fig[1, 2]; xlabel="mean temperature (K)", title="⟨T⟩(z)")
lines!(ax2, T̄, zc; linewidth=2.5); lines!(ax2, Tr, zr; linewidth=2.5, linestyle=:dash, color=:black)
ax3 = Axis(fig[1, 3]; xlabel="rms vertical velocity (m/s)", title="w′(z)")
lines!(ax3, sqrt.(ww), zc; linewidth=2.5); lines!(ax3, wr, zr; linewidth=2.5, linestyle=:dash, color=:black)
hideydecorations!(ax2; grid=false); hideydecorations!(ax3; grid=false)
save("figures/profiles_vs_reference.png", fig)
@info "saved profiles"

##### Calibration map
runs = [
    ("parity_long", "C 6e-3, dry", :dry), ("dry_C2e-2", "C 2e-2, dry", :dry), ("dry_C4e-2", "C 4e-2, dry", :dry),
    ("cloudy_tau72_C2e-2", "C 2e-2", :tau72), ("cloudy_tau72_C2.5e-2", "C 2.5e-2", :tau72), ("cloudy_tau72_C3e-2", "C 3e-2", :tau72), ("cloudy_tau72_C4e-2", "C 4e-2", :tau72),
    ("ladder128_cloudy_tau72_C2e-2", "C 2e-2, 128³", :tau72),
    ("cloudy_tau28_C2e-2", "C 2e-2", :tau28), ("cloudy_tau28_C2.5e-2", "C 2.5e-2", :tau28), ("cloudy_tau28_C3e-2", "C 3e-2", :tau28), ("cloudy_tau28_C3e-2_side1.5e-2", "C 3e-2 / sides 1.5e-2", :tau28),
    ("cloudy_tau28_loglaw1mm", "ℓ 1 mm, WENO9", :loglaw9), ("cloudy_tau28_loglaw0.5mm", "ℓ 0.5 mm, WENO9", :loglaw9),
    ("weno5_loglaw1mm", "ℓ 1 mm", :loglaw5), ("weno5_loglaw1.25mm", "ℓ 1.25 mm (default)", :loglaw5), ("weno5_loglaw1.5mm", "ℓ 1.5 mm", :loglaw5),
]
families = Dict(:dry => ("dry chamber, constant C", :gray55), :tau72 => ("cloudy, 72 s relaxation, constant C", :steelblue),
                :tau28 => ("cloudy, 28 s relaxation, constant C", :dodgerblue4), :loglaw9 => ("cloudy, 28 s, log law, WENO(9)", :darkorange3),
                :loglaw5 => ("cloudy, 28 s, log law, WENO(5)", :firebrick))
fig = Figure(size=(900, 620), fontsize=15)
ax = Axis(fig[1, 1]; xlabel="Lagrangian spread of 𝒮 along droplets (%)", ylabel="activated fraction after 60 s at a −1 % target (fluctuating)",
          title="Every chamber run against the reference replay")
seen = Set{Symbol}()
for (file, label, fam) in runs
    path = "runs/$(file)_activation.jld2"
    isfile(path) || continue
    f = load(path); results = f["results"]; series = f["series"]
    windows = sort([k for k in keys(results) if startswith(k, "window_") && haskey(series[k], "S")])
    isempty(windows) && continue                            # runs before the Lagrangian sampling existed
    m = findfirst(t -> isapprox(t, -0.01; atol=1e-6), results["targets"])
    act = mean(results[w]["fluctuating"][m] for w in windows)
    σ = mean(100 * std(series[w]["S"]) for w in windows)
    name, color = families[fam]
    scatter!(ax, [σ], [act]; color, markersize=13, label=(fam in seen ? nothing : name))
    push!(seen, fam)
    text!(ax, σ, act; text=label, fontsize=10, offset=(7, -4), align=(:left, :top), color)
end
ref = load("runs/reference_replay.jld2"); rr = ref["results"]; rs = ref["series"]
wr = sort([k for k in keys(rr) if startswith(k, "window_")])
m = findfirst(t -> isapprox(t, -0.01; atol=1e-6), rr["targets"])
ref_act = mean(rr[w]["fluctuating"][m] for w in wr); ref_σ = mean(100 * std(rs[w]["S"]) for w in wr)
vlines!(ax, [ref_σ]; color=:black, linestyle=:dash); hlines!(ax, [ref_act]; color=:black, linestyle=:dash)
scatter!(ax, [ref_σ], [ref_act]; color=:black, marker=:xcross, markersize=22, label="SAM replay")
axislegend(ax; position=:lt, labelsize=11)
save("figures/calibration_map.png", fig)
@info "saved calibration map" ref_σ ref_act
