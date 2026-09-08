# Animate the droplet population from `scripts/run_movie.jl`: the droplets near a vertical slice
# coloured by size and marked when activated, the evolving size spectrum, and the (𝒮, D) plane
# with the Köhler equilibrium curve, on which activation is the crossing of the critical point.
#
# Usage: julia --project analysis/make_particle_movie.jl [prefix] [output.mp4] [fps]
using Breeze, CairoMakie, JLD2, Printf, Statistics
using Breeze.LagrangianMicrophysics: equilibrium_supersaturation, critical_diameter, critical_supersaturation

prefix = length(ARGS) ≥ 1 ? ARGS[1] : "runs/movie"
output = length(ARGS) ≥ 2 ? ARGS[2] : "figures/droplets.mp4"
fps = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 20
stride = parse(Int, get(ENV, "MOVIE_STRIDE", "1"))
width = parse(Int, get(ENV, "MOVIE_WIDTH", "1300"))

file = jldopen("$(prefix)_droplets.jld2")
iters = sort(parse.(Int, filter(k -> k != "serialized", keys(file["timeseries/t"]))))
times = [file["timeseries/t/$i"] for i in iters]
frames = [file["timeseries/droplets/$i"] for i in iters]
close(file)
t₀ = times[1]
FT = Float64
constants = ThermodynamicConstants(FT)
Dᵈ, κ, T = FT(frames[1].Dᵈ[1]), FT(frames[1].κ[1]), FT(287.6)
Dᶜ = critical_diameter(Dᵈ, κ, T, constants)
𝒮ᶜ = critical_supersaturation(Dᵈ, κ, T, constants)
@info "droplet frames" length(frames) length(frames[1].x) 1e6Dᶜ 100𝒮ᶜ

# The Köhler equilibrium curve of this aerosol
Dcurve = 10 .^ range(log10(1.05 * Dᵈ), log10(20e-6), length=400)
𝒮curve = [100 * equilibrium_supersaturation(D, Dᵈ, κ, T, constants) for D in Dcurve]

diameters(p) = 1e6 .* sqrt.(max.(p.D², 0))
activated(p) = p.D² .≥ p.Dᶜ .^ 2
Dmax = maximum(maximum(diameters(p)) for p in frames)
𝒮all = reduce(vcat, [100 .* p.𝒮 for p in frames])
𝒮lim = (quantile(𝒮all, 0.002), quantile(𝒮all, 0.998))

n = Observable(1)
title = @lift @sprintf("%d droplets in the chamber   t = %5.1f s   activated %.1f %%",
                       length(frames[1].x), times[$n] - t₀, 100 * count(activated(frames[$n])) / length(frames[$n].x))

fig = Figure(size=(width, round(Int, 0.44 * width)), fontsize=round(Int, 15 * width / 1300))
Label(fig[0, 1:3], title; fontsize=round(Int, 19 * width / 1300), font=:bold, tellwidth=false)

# Left: the droplets near the slice, coloured by diameter
ax1 = Axis(fig[1, 1]; xlabel="x (m)", ylabel="z (m)", title="droplets within one cell of the slice", aspect=DataAspect())
near = @lift abs.(frames[$n].y) .< 2 / 32
positions = @lift Point2f.(frames[$n].x[$near], frames[$n].z[$near])
colors = @lift log10.(max.(diameters(frames[$n])[$near], 0.2))
sizes = @lift 2 .+ 6 .* (log10.(max.(diameters(frames[$n])[$near], 0.2)) .- log10(0.2)) ./ (log10(Dmax) - log10(0.2))
sc = scatter!(ax1, positions; color=colors, colormap=:turbo, colorrange=(log10(0.2), log10(Dmax)), markersize=sizes, strokewidth=0)
ticks = [0.3, 1, 3, 10, 30]
Colorbar(fig[2, 1], sc; vertical=false, flipaxis=false, label="diameter (μm)", height=12,
         ticks=(log10.(ticks), string.(ticks)))
xlims!(ax1, -1, 1); ylims!(ax1, 0, 1)

# Middle: the size spectrum
ax2 = Axis(fig[1, 2]; xlabel="diameter (μm)", ylabel="droplets", title="size spectrum", xscale=log10,
           limits=((0.15, 1.2Dmax), (0.7, 1.5 * length(frames[1].x))), yscale=log10)
edges = 10 .^ range(log10(0.15), log10(1.2Dmax), length=45)
centres = sqrt.(edges[1:end-1] .* edges[2:end])
counts = @lift begin
    h = zeros(length(centres))
    for D in diameters(frames[$n])
        k = searchsortedlast(edges, D)
        1 ≤ k ≤ length(h) && (h[k] += 1)
    end
    max.(h, 0.7)
end
barplot!(ax2, centres, counts; color=(:dodgerblue3, 0.8), gap=0.05)
vlines!(ax2, [1e6Dᶜ]; color=:firebrick, linestyle=:dash, label="critical diameter")
axislegend(ax2; position=:rt, labelsize=11)

# Right: the (𝒮, D) plane with the Köhler curve
ax3 = Axis(fig[1, 3]; xlabel="supersaturation seen by the droplet (%)", ylabel="diameter (μm)",
           title="Köhler plane", yscale=log10, limits=(𝒮lim, (0.15, 1.2Dmax)))
lines!(ax3, 𝒮curve, 1e6 .* Dcurve; color=:black, linewidth=2, label="equilibrium curve")
scatter!(ax3, [100𝒮ᶜ], [1e6Dᶜ]; color=:firebrick, markersize=12, marker=:diamond, label="critical point")
kohler = @lift Point2f.(100 .* frames[$n].𝒮, diameters(frames[$n]))
kcolor = @lift [a ? (:firebrick, 0.5) : (:steelblue, 0.35) for a in activated(frames[$n])]
scatter!(ax3, kohler; color=kcolor, markersize=4, strokewidth=0)
axislegend(ax3; position=:rb, labelsize=11)
vlines!(ax3, [0.0]; color=:gray70, linestyle=:dot)

colgap!(fig.layout, 1, 26); colgap!(fig.layout, 2, 26); rowgap!(fig.layout, 1, 6)
colsize!(fig.layout, 1, Relative(0.42))
record(fig, output, 1:stride:length(frames); framerate=fps) do i
    n[] = i
end
@info "saved" output
