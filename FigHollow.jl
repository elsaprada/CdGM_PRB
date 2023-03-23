#region ## LOAD CODE (MULTICORE)
include("code.jl")

using Distributed, ClusterManagers
const NPROCS = 32   # number of parallel processes per node
const NNODES = 6    # number of nodes
const NWORKERS = NNODES * NPROCS
if nprocs() <= NWORKERS
    addprocs_slurm(NWORKERS + 1 - nprocs(), retry_delays = Iterators.repeated(0.5),
    p = "esbirro", x  = "es[1]", ntasks_per_node = NPROCS,
    exeflags = "--project", enable_threaded_blas = false,  epilog = "./slurm.out/epilog.sh", job_file_loc = "slurm.out") 
end

# using Distributed
# const NWORKERS = 10   # number of parallel workers
# if nprocs() <= NWORKERS
#     addprocs(NWORKERS + 1 - nprocs())
# end

@everywhere begin
    using LinearAlgebra
    # using Pkg; Pkg.activate(".")
    BLAS.set_num_threads(1)
    include("code.jl")
end

#endregion

#region ## CONFIGURATION FIG 4 - COMMON

base_model = (;
    ξd = 70,  Rshell = 70, Rcore = 70, Raxis = 69, a0 = 5,
    Δ0 = 0.23, α=0, g = 0, Vexponent = 2)

base_params = (; τ = 0.0, τΓ = 30.0)

ωrng = range(-0.26, 0.26, 301) .+ 0.001im
nrng = range(0, 2.499, 400)

# ωrng = range(0.0, 0.26, 51) .+ 0.001im
# nrng = range(0, 2.499, 60)

function density_of_states()
    xs, ys = nrng, ωrng;
    mjdict = Dict([(@show Z; Z => dos(nω_leadS(Z), xs, ys)) for Z in mjs]);
    return xs, ys, mjdict
end

function plot((xs, ys, mjdict), filename;  zfactor = 1.0, zfactor´ = 1.0, xlabel = "Φ/Φ₀", labelbar = "LDOS [a.u]")
    plot_and_export(xs, ys, mjdict; filename, zfactor, zfactor´, xlabel, labelbar)
end

function bands(Z; n = 1, τΓ = 0.0, τ = 0.0, kmax = pi)
    hh = leadS(; Z, ω = 0.0000im, n, τΓ, τ)
    b = bandstructure(hh, cuboid((0, kmax), subticks = 81), method = LinearAlgebraPackage())
    vlplot(b, size = 400, ylims = (-40,40))
end

#endregion


#region ## PANEL 4(a)

model = (; base_model..., Raxis = 69, Vmin = -133.6, Vmax = -133.6)
params = (; base_params..., τΓ = 0.1)
mjs = -4:4

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τ = 0.0, τΓ = params.τΓ)

dataA = density_of_states()
plot(dataA, "Output/FigHollow/FigHollowa"; zfactor = 0.7, zfactor´ = 0.7)

# bands(0, kmax = 0.5)

# hh = leadS(; Z = 5, ω = 0+0.0001im, n = 1.0, τ = 0.0, τΓ = 0.0)
# vlplot(hh, maxthickness = 0.1)
#endregion

#region ## PANEL 4(b)

model = (; base_model..., Raxis = 69, Vmin = -133.6, Vmax = -133.6)
params = (; base_params..., τΓ = 0.8)
mjs = -4:4

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τ = 0.0, τΓ = params.τΓ)

dataB = density_of_states()
plot(dataB, "Output/FigHollow/FigHollowb"; zfactor = 0.5, zfactor´ = 0.7)

#endregion

#region ## PANEL 4(c)

model = (; base_model..., Raxis = 69, Vmin = -133.6, Vmax = -133.6)
params = (; base_params..., τΓ = 3)
mjs = -4:4

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τ = 0.0, τΓ = params.τΓ)

dataC = density_of_states()
plot(dataC, "Output/FigHollow/FigHollowc"; zfactor = 0.6, zfactor´ = 0.7)

#endregion

## Bands
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
bands(0; kmax = 0.4, n = 0, τΓ = 0)