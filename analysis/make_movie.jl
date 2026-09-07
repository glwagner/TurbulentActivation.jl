# Animate the chamber from the slices written by `scripts/run_movie.jl`: the speed and the
# supersaturation on a vertical slice through the middle of the box, with the droplets that lie
# within half a cell of that slice drawn on top, activated ones filled.
#
# Usage: julia --project analysis/make_movie.jl [prefix] [output.mp4 or .gif] [fps]
using Breeze, CairoMakie, JLD2, Printf, Statistics
using Oceananigans
using Oceananigans.Fields: interior

prefix = length(ARGS) ≥ 1 ? ARGS[1] : "runs/movie"
output = length(ARGS) ≥ 2 ? ARGS[2] : "figures/chamber.mp4"
fps = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 20

speed = FieldTimeSeries("$(prefix)_slices.jld2", "speed")
𝒮 = FieldTimeSeries("$(prefix)_slices.jld2", "𝒮")
times = speed.times
grid = speed.grid
x = xnodes(speed[1]) ; z = znodes(speed[1])
@info "frames" length(times) extrema(times) size(interior(speed[1]))

droplets = nothing
dfile = "$(prefix)_droplets.jld2"
if isfile(dfile)
    droplets = jldopen(dfile)
    @info "droplet frames" length(keys(droplets["timeseries/t"]))
end

smax = maximum(maximum(abs, interior(speed[n])) for n in 1:length(times))
slim = (0, round(smax, digits=2))
𝒮lim = let a = maximum(maximum(abs, interior(𝒮[n])) for n in 1:length(times))
    (-100a, 100a)
end
@info "ranges" slim 𝒮lim

n = Observable(1)
speedₙ = @lift interior(speed[$n], :, 1, :)
𝒮ₙ = @lift 100 .* interior(𝒮[$n], :, 1, :)
title = @lift @sprintf("Pi Chamber, vertical slice   t = %5.1f s", times[$n] - times[1])

fig = Figure(size=(1150, 460), fontsize=16)
Label(fig[0, 1:2], title; fontsize=19, font=:bold)
ax1 = Axis(fig[1, 1]; xlabel="x (m)", ylabel="z (m)", title="speed (m/s)", aspect=DataAspect())
hm1 = heatmap!(ax1, x, z, speedₙ; colormap=:speed, colorrange=slim)
Colorbar(fig[1, 1, Right()], hm1; width=12)
ax2 = Axis(fig[1, 2]; xlabel="x (m)", title="supersaturation (%)", aspect=DataAspect())
hm2 = heatmap!(ax2, x, z, 𝒮ₙ; colormap=:balance, colorrange=𝒮lim)
Colorbar(fig[1, 2, Right()], hm2; width=12)
hideydecorations!(ax2; grid=false)

if !isnothing(droplets)
    iters = sort(parse.(Int, filter(k -> k != "serialized", keys(droplets["timeseries/t"]))))
    Δy = 2 / size(grid, 2)                       # the slice is one cell thick
    frame_positions = map(iters) do i
        p = droplets["timeseries/droplets/$i"]
        near = abs.(p.y) .< Δy
        (Point2f.(p.x[near], p.z[near]), Float32.(1e6 .* sqrt.(p.D²[near])))
    end
    positions = @lift frame_positions[min($n, length(frame_positions))][1]
    sizes = @lift 1.5 .+ 2.5 .* log10.(max.(frame_positions[min($n, length(frame_positions))][2], 0.1f0) .+ 1)
    for ax in (ax1, ax2)
        scatter!(ax, positions; markersize=sizes, color=(:white, 0.55), strokecolor=(:black, 0.45), strokewidth=0.3)
    end
end

record(fig, output, 1:length(times); framerate=fps) do i
    n[] = i
end
@info "saved" output
