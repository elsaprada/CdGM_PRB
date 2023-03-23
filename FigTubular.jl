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

@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    include("code.jl")
end

#endregion

#region ## CONFIGURATION FIG 5 - COMMON

base_model = (;
    ξd = 70, Rshell = 70, Rcore = 70, Raxis = 0, a0 = 5,
    Δ0 = 0.23, α=0, g = 0, Vexponent = 2)

base_params = (; τΓ = 10)

ωrng = range(-0.26, 0.26, 300) .+ 0.001im
nrng = range(0, 2.499, 400)

# ωrng = range(0, 0.26, 41) .+ 0.001im
# nrng = range(0, 2.499, 40)

function density_of_states()
    xs, ys = nrng, ωrng;
    mjdict = Dict([(@show Z; Z => dos(nω_leadS(Z), xs, ys)) for Z in mjs]);
    return xs, ys, mjdict
end

function plot((xs, ys, mjdict), filename;  zfactor = 1.0, zfactor´ = 1.0, xlabel = "Φ/Φ₀", labelbar = "LDOS [a.u]")
    plot_and_export(xs, ys, mjdict; filename, zfactor, zfactor´, xlabel, labelbar)
end

function bands(Z; nev = 40, n = 1, τΓ = 0.0, kmax = pi)
    hh = leadS(; Z, ω = 0, n, τΓ)
    b = bandstructure(hh, cuboid((0, kmax), subticks = 81), method = ArpackPackage(; nev, sigma = 0.1im))
    vlplot(b, size = 400, ylims = (-40,40))
end

#endregion

#region ## PANEL 5(b)

model = (; base_model..., Raxis = 70, Vmin = -145, Vmax = -145)
params = (; base_params..., τΓ = 3)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)
# bands(6; nev = 2)

dataB = density_of_states()
# plot(dataB, "Output/FigTubular/FigTubularb"; zfactor = 0.7, zfactor´ = 0.9)
#endregion

#region ## PANEL 5(c)

model = (; base_model..., Raxis = 60, Vmin = -52, Vmax = -52)
params = (; base_params..., τΓ = 12)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

# bands(6)

dataC = density_of_states()
# plot(dataC, "Output/FigTubular/FigTubularc"; zfactor = 0.7, zfactor´ = 0.6)


#endregion

#region ## PANEL 5(d)

model = (; base_model..., Raxis = 50, Vmin = -34, Vmax = -34)
params = (; base_params..., τΓ = 40)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

# bands(6)

dataD = density_of_states()
# plot(dataD, "Output/FigTubular/FigTubulard"; zfactor = 0.5, zfactor´ = 0.6)

#endregion

#region ## PANEL 5(e)

model = (; base_model..., Raxis = 40, Vmin = -29, Vmax = -29)
params = (; base_params..., τΓ = 90)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

# bands(6)

dataE = density_of_states()
# plot(dataE, "Output/FigTubular/FigTubulare"; zfactor = 0.5, zfactor´ = 0.6)

#endregion

#region ## PANEL 5(f)

model = (; base_model..., Raxis = 30, Vmin = -28, Vmax = -28)
params = (; base_params..., τΓ = 170)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

# bands(7)

dataF = density_of_states()
# plot(dataF, "Output/FigTubular/FigTubularf"; zfactor = 0.8, zfactor´ = 0.6)

#endregion

#region ## PANEL 5(g)

model = (; base_model..., Raxis = 20, Vmin = -28, Vmax = -28)
params = (; base_params..., τΓ = 350)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

# bands(6)

dataG = density_of_states()
# plot(dataG, "Output/FigTubular/FigTubularg"; zfactor = 0.8, zfactor´ = 0.6)


#endregion

#region ## PANEL 5(h)

model = (; base_model..., Raxis = 10, Vmin = -28, Vmax = -28)
params = (; base_params..., τΓ = 350)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

# bands(6)

dataH = density_of_states()
# plot(dataH, "Output/FigTubular/FigTubularh"; zfactor = 0.5, zfactor´ = 0.6)

#endregion

#region ## PANEL 5(i)

model = (; base_model..., Raxis = 0, Vmin = -28, Vmax = -28)
params = (; base_params..., τΓ = 350)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

dataI = density_of_states()
# plot(dataI, "Output/FigTubular/FigTubulari"; zfactor = 0.5, zfactor´ = 0.6)

#endregion


#region ## Plot panels

plot(dataB, "Output/FigTubular/FigTubularb"; zfactor = 0.5, zfactor´ = 0.6);
plot(dataC, "Output/FigTubular/FigTubularc"; zfactor = 0.5, zfactor´ = 0.6);
plot(dataD, "Output/FigTubular/FigTubulard"; zfactor = 0.5, zfactor´ = 0.6);
plot(dataE, "Output/FigTubular/FigTubulare"; zfactor = 0.5, zfactor´ = 0.6);
plot(dataF, "Output/FigTubular/FigTubularf"; zfactor = 0.5, zfactor´ = 0.6);
plot(dataG, "Output/FigTubular/FigTubularg"; zfactor = 0.5, zfactor´ = 0.6);
plot(dataH, "Output/FigTubular/FigTubularh"; zfactor = 0.5, zfactor´ = 0.6);
plot(dataI, "Output/FigTubular/FigTubulari"; zfactor = 0.5, zfactor´ = 0.6);

#endregion
