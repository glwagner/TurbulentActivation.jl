# Vertical profiles of the chamber at several vertical resolutions against the reference SAM
# fields: the wall layers are one cell thick at the reference's own spacing, so the mean state
# depends on Δz even where the fluctuation statistics do not.
#
# Usage: julia --project analysis/resolution_profiles.jl out.png label=prefix [label=prefix ...]
using JLD2, CairoMakie, Statistics, Printf
include(joinpath(@__DIR__, "..", "reference", "read_bin3d.jl"))
using .SAMBin3D

output = ARGS[1]
runs = [(label = rsplit(a, "="; limit=2)[1], prefix = rsplit(a, "="; limit=2)[2]) for a in ARGS[2:end]]

function profiles(prefix)
    d = jldopen("runs/$(prefix)_profiles.jld2")
    its = sort(parse.(Int, filter(k -> k != "serialized", keys(d["timeseries/T"]))))
    window = its[max(1, end - 36):end]
    mean_profile(k) = mean([vec(d["timeseries/$k/$i"]) for i in window])
    T, ℋ, ww = mean_profile("T"), mean_profile("ℋ"), mean_profile("ww")
    close(d)
    keep = T .> 0
    T, ℋ, ww = T[keep], ℋ[keep], ww[keep]
    Nz = length(T)
    return (z = ((1:Nz) .- 0.5) ./ Nz, T = T, 𝒮 = 100 .* (ℋ .- 1), w = sqrt.(ww))
end

# Reference
dir = joinpath(homedir(), "anderson_reference", "LES", "LES")
esat(T) = 610.94 * exp(17.625 * (T - 273.15) / (T - 30.11))
Tr, Sr, wr = zeros(32), zeros(32), zeros(32)
nsnap = 0
for step in 90000:15000:180000
    s = read_bin3d(joinpath(dir, @sprintf("PiChamber_huji_19K_trj_32_%010d.bin3D", step)))
    q = s.QV ./ 1000
    S = @. 1e5 * q / (0.622 + q) / esat(s.TABS) - 1
    for k in 1:32
        Tr[k] += mean(view(s.TABS, :, :, k)); Sr[k] += mean(view(S, :, :, k)); wr[k] += mean(view(s.W, :, :, k) .^ 2)
    end
    global nsnap += 1
end
reference = (z = ((1:32) .- 0.5) ./ 32, T = Tr ./ nsnap, 𝒮 = 100 .* Sr ./ nsnap, w = sqrt.(wr ./ nsnap))

fig = Figure(size=(1050, 460), fontsize=15)
colors = [:dodgerblue3, :darkorange2, :seagreen4, :purple3]
axes = []
for (i, (field, name, unit)) in enumerate(((:𝒮, "⟨𝒮⟩", "%"), (:T, "⟨T⟩", "K"), (:w, "w′", "m/s")))
    ax = Axis(fig[1, i]; xlabel="$name ($unit)", ylabel=i == 1 ? "z / H" : "", title=name * "(z)")
    for (j, run) in enumerate(runs)
        p = profiles(run.prefix)
        lines!(ax, getproperty(p, field), p.z; linewidth=2.5, color=colors[mod1(j, end)], label=run.label)
    end
    lines!(ax, getproperty(reference, field), reference.z; linewidth=2.5, linestyle=:dash, color=:black, label="SAM reference")
    field === :𝒮 && vlines!(ax, [0.0]; color=:gray60, linestyle=:dot)
    i > 1 && hideydecorations!(ax; grid=false)
    push!(axes, ax)
end
axislegend(axes[1]; position=:rb, labelsize=11)
save(output, fig)
@info "saved" output
for run in runs
    p = profiles(run.prefix)
    @printf("%-18s ⟨T⟩ = %.3f K (reference %.3f)   ⟨𝒮⟩ = %+.3f %% (%+.3f)   near-floor ΔT = %.2f K (%.2f)\n",
            run.label, mean(p.T), mean(reference.T), mean(p.𝒮), mean(reference.𝒮),
            p.T[1] - mean(p.T[3:end-2]), reference.T[1] - mean(reference.T[3:end-2]))
end
