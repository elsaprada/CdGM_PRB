#region ## LOAD CODE (MULTICORE)
include("code.jl")

using Distributed, ClusterManagers
const NPROCS = 32   # number of parallel processes per node
const NNODES = 5    # number of nodes
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

#region ## FUNCTIONS

function density_of_states()
    xs, ys = nrng, ωrng;
    mjdict = Dict([(@show Z; Z => dos(nω_leadS(Z), xs, ys)) for Z in mjs]);
    return xs, ys, mjdict
end

function local_density_of_states(x)
    ys = ωrng;
    mjdict = Dict([(@show Z; Z => dos_r(nω_leadS(Z), x, ys)) for Z in mjs]);
    xs = last.(Quantica.allsitepositions(leadS().lattice))
    return xs, ys, mjdict
end

function gamma()
    xs, ys = τrng, ωrng;
    mjdict = Dict([(@show Z; Z => dos(Γω_leadS(Z), xs, ys)) for Z in mjs]);
    return xs, ys, mjdict
end

function dosL()
    xs, ys, La0 = nrng, ωrng, model.L/model.a0;
    mjdict = Dict([(@show Z; Z => dos(nω_leadS(Z), xs, ys, La0, leadS)) for Z in mjs]);
    return xs, ys, mjdict
end

function conductance()
    xs, ys = nrng, ωrng;
    mjdict = Dict([(@show Z; Z => conductance(Lsys, ΓL, xs, ys; params..., Z)) for Z in mjs]);
    return xs, ys, mjdict
end

function plot((xs, ys, mjdict), filename;  zfactor = 1.0, zfactor´ = 1.0, xlabel = "Φ/Φ₀", labelbar = "LDOS [a.u]", kw...)
    plot_and_export(xs, ys, mjdict; filename, zfactor, zfactor´, xlabel, labelbar, kw...)
end

function bands(Z; nev = 40, n = 1, τ = 0.0, kmax = pi)
    hh = leadS(; Z, ω = 0, n, τ)
    b = bandstructure(hh, cuboid((0, kmax), subticks = 81), method = ArpackPackage(; nev, sigma = 0.1im))
    vlplot(b, size = 400, ylims = (-40,40))
end

#endregion

#region ## CONFIGURATION FIG SOC - COMMON (like FIG SOLID)
# ξ=70nm, Rshell = 80nm, dshell = 10nm
base_model = (;
    ξd = 70, Rshell = 80, Rcore = 70, Raxis = 0, a0 = 5, Vexponent = 2,
    Vmin = -150,
    Δ0 = 0.23,
    α = 0, α0 = 0, g = 0,
    L = 85, δL = 80, δVL = 30, δR = 0, δVR = 0)

base_params = (; n = 1, τΓ = 0)

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ);

τrng = range(0, 2, 30)

ωrng = range(0, 0.26, 201) .+ 0.00005im
nrng = range(0, 2.499, 400)

# ωrng = range(-0.26, 0.26, 101) .+ 0.001im
# nrng = range(0, 2.499, 150)

# ωrng = range(-0.26, 0.26, 81) .+ 0.00005im
# nrng = range(0, 2.499, 50)

#endregion

#region ## Fig SOC - top
modelSOC = (; base_model..., Vmin = -35, Vmax = 0, a0 = 5,
         Vmin_lead = -35, Vmax_lead = 0,
         δVL = 80, δL = 50, δSL = 50, L = 50,
         α = 1.19)
modelNOSOC = (; modelSOC..., α = 0)

params = (; base_params..., n = 0, τΓ = 90)

mjs = -7:7



## LDOS SOC
model = modelSOC
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosn = density_of_states()
plot(ldosn, "Output/FigSOC/FigSOC_top_ldos"; zfactor = 0.6)

## LDOS NOSOC
model = modelNOSOC
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosn0 = density_of_states()
plot(ldosn0, "Output/FigSOC/FigNOSOC_top_ldos"; zfactor = 0.2)

#endregion

#region ## Fig SOC - no-top
modelSOC = (; base_model..., Vmin = -40, Vmax = 0, a0 = 5,
         Vmin_lead = -40, Vmax_lead = 0,
         δVL = 80, δL = 50, δSL = 50, L = 50,
         α = 1.19)
modelNOSOC = (; modelSOC..., α = 0)

params = (; base_params..., n = 0, τΓ = 90)

mjs = -7:7

ωrng = range(0, 0.26, 201) .+ 0.00001im
nrng = range(0, 2.499, 400)

## LDOS SOC
model = modelSOC
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosn = density_of_states()
plot(ldosn, "Output/FigSOC/FigSOC_notop_ldos"; zfactor = 0.2)

## LDOS NOSOC
model = modelNOSOC
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosn0 = density_of_states()
plot(ldosn0, "Output/FigSOC/FigNOSOC_notop_ldos"; zfactor = 0.2)

#endregion

## Lattice
hh = leadS(; Z = 8, ω = 0+0.0001im, n = 0, τ = 1)
vlplot(hh, maxthickness = 0.06)

## Lattice
hh = Lsys(; Z = 8, ω = 0.0+0.0001im, n = 1, τ = 0.0)
vlplot(hh, maxthickness = 0.01)

## Bands
model = modelSOC
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
bands(0; nev = 12, τΓ = 0, n = 1)

## Load and plot
dict = load("Output/FigSOC/FigNOSOC_notop_ldos.jld2")
xs, ys, zs = dict["xs"], dict["ys"], dict["mjdict"]
params, model = dict["params"], dict["model"]
plot_and_export(xs, ys, zs; filename =  "Output/FigSOC/FigSOC_test", labelbar = "dI/dV [G₀]", zfactor = 0.9, xlabel = "Φ/Φ₀", pages = false, layers = false)