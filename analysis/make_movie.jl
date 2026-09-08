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
stride = parse(Int, get(ENV, "MOVIE_STRIDE", "1"))       # use every stride-th frame
width = parse(Int, get(ENV, "MOVIE_WIDTH", "1250"))      # figure width in points

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

# Robust colour ranges: the wall cells reach several times the interior values, so the ranges are
# set by a high quantile of all frames and the extremes are clipped
allspeed = reduce(vcat, vec(interior(speed[n])) for n in 1:length(times))
all𝒮 = reduce(vcat, vec(interior(𝒮[n])) for n in 1:length(times))
smax = round(quantile(allspeed, 0.995), digits=2)
𝒮max = round(100 * quantile(abs.(all𝒮), 0.99), digits=1)
@info "colour ranges" smax 𝒮max extrema(allspeed) 100 .* extrema(all𝒮)

n = Observable(1)
speedₙ = @lift interior(speed[$n], :, 1, :)
𝒮ₙ = @lift 100 .* interior(𝒮[$n], :, 1, :)
title = @lift @sprintf("Pi Chamber, vertical slice through the middle of the box   t = %5.1f s", times[$n] - times[1])

fig = Figure(size=(width, round(Int, 0.416 * width)), fontsize=round(Int, 16 * width / 1250))
Label(fig[0, 1:2], title; fontsize=19, font=:bold, tellwidth=false)
ax1 = Axis(fig[1, 1]; xlabel="x (m)", ylabel="z (m)", title="speed", aspect=DataAspect())
hm1 = heatmap!(ax1, x, z, speedₙ; colormap=:tempo, colorrange=(0, smax), highclip=:black)
Colorbar(fig[2, 1], hm1; vertical=false, flipaxis=false, label="m s⁻¹", height=12)
ax2 = Axis(fig[1, 2]; xlabel="x (m)", title="supersaturation", aspect=DataAspect())
hm2 = heatmap!(ax2, x, z, 𝒮ₙ; colormap=:balance, colorrange=(-𝒮max, 𝒮max), lowclip=:midnightblue, highclip=:darkred)
Colorbar(fig[2, 2], hm2; vertical=false, flipaxis=false, label="%", height=12)
hideydecorations!(ax2; grid=false)
rowgap!(fig.layout, 1, 6)
colgap!(fig.layout, 1, 24)

record(fig, output, 1:stride:length(times); framerate=fps) do i
    n[] = i
end
@info "saved" output
