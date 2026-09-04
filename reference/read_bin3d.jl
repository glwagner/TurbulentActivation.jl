# Read a SAM `bin3D` snapshot (Fortran sequential unformatted, see `bin3D2nc.f`) into a
# NamedTuple of fields on the (x, y, z) grid. Each record is framed by 4-byte length markers;
# the field blocks of this data set are single precision but twice the subdomain size,
# with the trailing half unused (the writer dumped an oversized buffer; 9.4 MB per file).
#
#   snapshot = read_bin3d("PiChamber_huji_19K_trj_32_0000120000.bin3D")
#   snapshot.TABS, snapshot.QV, snapshot.z, snapshot.p, snapshot.time
module SAMBin3D

export read_bin3d

function record(io)
    n = read(io, Int32)
    bytes = read(io, n)
    n′ = read(io, Int32)
    n == n′ || error("corrupt Fortran record: $n ≠ $n′")
    return bytes
end

reinterpret_record(io, T) = reinterpret(T, record(io))

function read_bin3d(filename)
    open(filename, "r") do io
        header = reinterpret_record(io, Int32)
        nx, ny, nz, nsubs, nsubsx, nsubsy, nfields = Int.(header[1:7])
        z = [reinterpret_record(io, Float64)[1] for _ in 1:nz]
        p = [reinterpret_record(io, Float64)[1] for _ in 1:nz]
        dx = reinterpret_record(io, Float64)[1]
        dy = reinterpret_record(io, Float64)[1]
        time = reinterpret_record(io, Float64)[1]
        nx_gl, ny_gl = nx * nsubsx, ny * nsubsy
        x = dx .* (0:nx_gl-1)
        y = dy .* (0:ny_gl-1)
        fields = Dict{Symbol, Array{Float32, 3}}()
        names = String[]
        for _ in 1:nfields
            label = record(io)                       # name(8) blank(1) long_name(80) blank(1) units(10)
            name = strip(String(label[1:8]))
            field = Array{Float32}(undef, nx_gl, ny_gl, nz)
            for n in 0:nsubs-1
                values = reinterpret_record(io, Float32)
                sub = reshape(values[1:nx*ny*nz], nx, ny, nz)
                i0 = nx * (n % nsubsx)
                j0 = ny * (n ÷ nsubsx)
                field[i0+1:i0+nx, j0+1:j0+ny, :] .= sub
            end
            fields[Symbol(name)] = field
            push!(names, name)
        end
        return (; x, y, z, p, time, dx, dy, names, fields...)
    end
end

end # module
