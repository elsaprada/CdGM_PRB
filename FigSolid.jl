#region ## LOAD CODE (MULTICORE)
include("code.jl")

using Distributed, ClusterManagers
const NPROCS = 32   # number of parallel processes per node
const NNODES = 3    # number of nodes
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
    mjdict = Dict([(@show Z; Z => conductance(Lsys, ΓL, xs, ys; params..., Z, slead = true)) for Z in mjs]);
    return xs, ys, mjdict
end

function plot((xs, ys, mjdict), filename;  zfactor = 1.0, zfactor´ = 1.0, xlabel = "Φ/Φ₀", labelbar = "LDOS [a.u]", kw...)
    plot_and_export(xs, ys, mjdict; filename, zfactor, zfactor´, xlabel, labelbar, kw...)
end

function bands(Z; nev = 40, kmax = pi, subticks = 81, params...)
    hh = leadS(; Z, params...)
    b = bandstructure(hh, cuboid((0, kmax); subticks), method = ArpackPackage(; nev, sigma = 0.1im))
    vlplot(b, size = 400, ylims = (-40,40))
end

#endregion

#region ## CONFIGURATION FIG 7 - COMMON
# ξ=70nm, Rshell = 80nm, dshell = 10nm
base_model = (;
    ξd = 70, Rshell = 80, Rcore = 70, Raxis = 0, a0 = 5, Vexponent = 2,
    Vmin = -150,
    Δ0 = 0.23,
    α = 0, α0 = 0, g = 0,
    L = 85, δL = 80, δVL = 30, δR = 0, δVR = 0)

base_params = (; n = 1, τΓ = 0)

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ);
Γω_leadS(Z) = (τΓ, ω) -> leadS(; Z = Z, ω = ω, n = params.n, τΓ)

τrng = range(0, 2, 30)

ωrng = range(-0.26, 0.26, 300) .+ 0.001im
nrng = range(0, 2.499, 400)

# ωrng = range(-0.26, 0.26, 101) .+ 0.001im
# nrng = range(0, 2.499, 150)

# ωrng = range(-0.26, 0.26, 80) .+ 0.001im
# nrng = range(0, 2.499, 50)

#endregion

#region ## Column A
colname = "A"
model0 = (; base_model..., Vmin = -40, Vmax = 0, a0 = 5,
         Vmin_lead = -40, Vmax_lead = 0)
model1 = (; model0..., δVL = 60, δL = 50, δSL = 50, L = 50)
model3 = (; model0..., δVL = 25, δL = 150, δSL = 150, L = 150)

params = (; base_params..., n = 0, τΓ = 90)

mjs = -4:4

## LDOS
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosA = density_of_states()
#plot(ldosA, "Output/FigSolid/FigSolid_$(colname)_ldos"; zfactor = 0.5)

## Cond short
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
condAs = conductance()
#plot(condAs, "Output/FigSolid/FigSolid_$(colname)_cond_short", labelbar = "dI/dV [G₀]", zfactor = 0.5)

## Cond long
model = model3
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
condAl = conductance()
#plot(condAl, "Output/FigSolid/Test/FigSolid_$(colname)_cond_long", labelbar = "dI/dV [G₀]", zfactor = 0.5)

## Transparency
mjs = 0
gm = gamma()
# plot(gm, "Output/FigSolid/FigSolid_$(colname)_gamma"; zfactor = 0.9, xlabel = "τ")


#endregion

#region ## Column B
colname = "B"
model0 = (; base_model..., Vmin = -70, Vmax = 0, a0 = 5,
         Vmin_lead = -70, Vmax_lead = 0)
model1 = (; model0..., δVL = 110, δL = 50, δSL = 50, L = 50)
model3 = (; model0..., δVL = 50, δL = 150, δSL = 150, L = 150)

params = (; base_params..., n = 0, τΓ = 30)

mjs = -7:7

## LDOS
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosB = density_of_states()
# plot(ldosB, "Output/FigSolid/FigSolid_$(colname)_ldos"; zfactor = 0.5)

## Cond short
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
condBs = conductance()
# plot(condBs, "Output/FigSolid/FigSolid_$(colname)_cond_short", labelbar = "dI/dV [G₀]", zfactor = 0.4)

## Cond long
model = model3
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
condBl = conductance()
# plot(condBl, "Output/FigSolid/FigSolid_$(colname)_cond_long", labelbar = "dI/dV [G₀]", zfactor = 0.5)

## Transparency
mjs = 0
gm = gamma()
plot(gm, "Output/FigSolid/FigSolid_$(colname)_gamma"; zfactor = 0.9, xlabel = "τ")


## Extra
model = (; model1..., Vmax = -30)
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
condBs = conductance()
plot(condBs, "Output/FigSolid/FigSolid_$(colname)_cond_short_extra", labelbar = "dI/dV [G₀]", zfactor = 0.4)

#endregion

#region ## Column C
colname = "C"
model0 = (; base_model..., Vmin = -140, Vmax = 0, a0 = 5,
         Vmin_lead = -140, Vmax_lead = 0)
model1 = (; model0..., δVL = 170, δL = 50, δSL = 50, L = 50)
model3 = (; model0..., δVL = 110, δL = 150, δSL = 150, L = 150)

Lsys, iC, iL, ΓL, leadS = build_cyl(; model1...);

params = (; base_params..., n = 0, τΓ = 20)

mjs = -12:12
# mjs = -4:4

## LDOS
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosC = density_of_states()
# plot(ldosC, "Output/FigSolid/FigSolid_$(colname)_ldos"; zfactor = 0.5)

## Cond short
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
condCs = conductance()
# plot(condCs, "Output/FigSolid/FigSolid_$(colname)_cond_short", labelbar = "dI/dV [G₀]", zfactor = 0.5)

## Cond long
model = model3
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
condCl = conductance()
# plot(condCl, "Output/FigSolid/FigSolid_$(colname)_cond_long", labelbar = "dI/dV [G₀]", zfactor = 0.5)

## Transparency
mjs = 0
gm = gamma()
plot(gm, "Output/FigSolid/FigSolid_$(colname)_gamma"; zfactor = 0.9, xlabel = "τ")

#endregion

## Plots

plot(ldosA, "Output/FigSolid/FigSolid_A_ldos"; zfactor = 0.5);
plot(condAs, "Output/FigSolid/FigSolid_A_cond_short", labelbar = "dI/dV [G₀]", zfactor = 0.5);
plot(condAl, "Output/FigSolid/FigSolid_A_cond_long", labelbar = "dI/dV [G₀]", zfactor = 0.5);

plot(ldosB, "Output/FigSolid/FigSolid_B_ldos"; zfactor = 0.5);
plot(condBs, "Output/FigSolid/FigSolid_B_cond_short", labelbar = "dI/dV [G₀]", zfactor = 0.5);
plot(condBl, "Output/FigSolid/FigSolid_B_cond_long", labelbar = "dI/dV [G₀]", zfactor = 0.5);

plot(ldosC, "Output/FigSolid/FigSolid_C_ldos"; zfactor = 0.5);
plot(condCs, "Output/FigSolid/FigSolid_C_cond_short", labelbar = "dI/dV [G₀]", zfactor = 0.5);
plot(condCl, "Output/FigSolid/FigSolid_C_cond_long", labelbar = "dI/dV [G₀]", zfactor = 0.5);

## Lattice
hh = leadS(; Z = 8, ω = 0+0.0001im, n = 0, τ = 1)
vlplot(hh, maxthickness = 0.06)

## Lattice
hh = Lsys(; Z = 8, ω = 0.0+0.0001im, n = 1, τ = 0.0)
vlplot(hh, maxthickness = 0.01)

## Bands
model = (; base_model..., Rcore = 300, Rshell = 300, Vmin = -1000, Vmax = -1000, Vs = -1000, a0 = 1, as = 1, m0 = 0.026, ms = 0.026, Δ = 100)
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
bands(0; ω = 0+0.000001im, nev = 22, τ = 1, n = 1, subticks = 181, kmax = 1.5, Δcore = 100)

## Load and plot
dict = load("Output/FigSolid/FigSolid_CZ_cond_long.jld2")
xs, ys, zs = dict["xs"], dict["ys"], dict["mjdict"]
params, model = dict["params"], dict["model"]
plot_and_export(xs, ys, zs; filename =  "Output/FigSolid/FigSolid_CZ_cond_long", labelbar = "dI/dV [G₀]", zfactor = 0.2, xlabel = "Φ/Φ₀", pages = false, layers = false)