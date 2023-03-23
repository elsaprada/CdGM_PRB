#region ## LOAD CODE (MULTICORE)
include("code.jl")

# using Distributed, ClusterManagers
# const NPROCS = 32   # number of parallel processes per node
# const NNODES = 5    # number of nodes
# const NWORKERS = NNODES * NPROCS
# if nprocs() <= NWORKERS
#     addprocs_slurm(NWORKERS + 1 - nprocs(), retry_delays = Iterators.repeated(0.5),
#     p = "esbirro", x  = "es[1]", ntasks_per_node = NPROCS,
#     exeflags = "--project", enable_threaded_blas = false,  epilog = "./slurm.out/epilog.sh", job_file_loc = "slurm.out") 
# end

using Distributed
const NWORKERS = 10   # number of parallel workers
if nprocs() <= NWORKERS
    addprocs(NWORKERS + 1 - nprocs())
end

@everywhere begin
    using LinearAlgebra
    using Pkg; Pkg.activate(".")
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

function bands(Z; nev = 40, n = 1, τΓ = 0.0, kmax = pi)
    hh = leadS(; Z, ω = 0, n, τΓ)
    b = bandstructure(hh, cuboid((0, kmax), subticks = 81), method = ArpackPackage(; nev, sigma = 0.1im))
    vlplot(b, size = 400, ylims = (-40,40))
end

#endregion

#region ## CONFIGURATION FIG 8 - COMMON
# ξ=70nm, Rshell = 80nm, dshell = 10nm
base_model = (;
    ξd = 70, Rshell = 80, Rcore = 70, Raxis = 0, a0 = 2, Vexponent = 2,
    Vmin = -150,
    Δ0 = 0.23,
    α = 0, α0 = 0, g = 0,
    L = 85, δL = 80, δVL = 30, δR = 0, δVR = 0)

base_params = (; n = 1, τΓ = 0)

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ);
Γω_leadS(Z) = (τ, ω) -> leadS(; Z = Z, ω = ω, n = params.n, τΓ)

τrng = range(0, 2, 30)

ωrng = range(0, 0.26, 151) .+ 0.001im
nrng = range(0.000001, 2.499, 301)
# nrng = range(0.0, 2.499, 151)

ωrng = range(0, 0.26, 41) .+ 0.001im
# nrng = range(0.0000, 2.499, 141)
nrng = range(0.00000, 2.499, 41)


#endregion

#region ## Fig 8ab - Similar to 7c, but with different Vmin and a0 (same # of occupied modes)
model0 = (; base_model..., Vmin = -70, Vmax = 0, a0 = 2,
         Vmin_lead = -70, Vmax_lead = 0)
model1 = (; model0..., δVL = 80, δL = 50, δSL = 50, L = 50)

params = (; base_params..., n = 0, τΓ = 600)

mjs = -7:7

## LDOS
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosn = density_of_states()
plot(ldosn, "Output/FigEffective/FigEffective_a_ldos_test"; zfactor = 0.9)

## LDOS vs r
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosr = local_density_of_states(1)
plot(ldosr, "Output/FigEffective/FigEffective_b_ldos"; zfactor = 0.3, layers = false)

#endregion

#region ## Fig 8c - Effective hollow
# model0 = (; base_model..., Rav = 59,
#          Vmin = -861.5, Vmax = -861.5, a0 = 2,
#          Vmin_lead = -861.5, Vmax_lead = -861.5)
# model1 = (; model0..., δVL = 80, δL = 50, δSL = 50, L = 50)

model0 = (; base_model..., Rav = 59,
         Vmin = -861.5+18, Vmax = -861.5+18, a0 = 2,
         Vmin_lead = -861.5+18, Vmax_lead = -861.5+18)
model1 = (; model0..., δVL = 80, δL = 50, δSL = 50, L = 50)

model0 = (; base_model..., Rav = 56,
         Vmin = -861.5+20, Vmax = -861.5+20, a0 = 2,
         Vmin_lead = -861.5+20, Vmax_lead = -861.5+20)
model1 = (; model0..., δVL = 80, δL = 50, δSL = 50, L = 50)


params = (; base_params..., n = 0, τΓ = 4.3)

mjs = -7:7

## LDOS
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosn = density_of_states()
plot(ldosn, "Output/FigEffective/FigEffective_c_ldos"; zfactor = 0.75)

## LDOS vs r
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
ldosr = local_density_of_states(1)
plot(ldosr, "Output/FigEffective/FigEffective_d_ldos"; zfactor = 0.9, layers = false)

#endregion

#endregion

## Lattice
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
hh = leadS(; Z = 0, ω = 0.3+0.0001im, n = 0, τΓ = 10)
vlplot(hh, maxthickness = 0.06)

## Bands
model = model1
Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);
bands(4; kmax = 0.4, nev = 2, n = 0, τΓ = 0)

## Load and plot
dict = load("Output/FigEffective/FigEffectivea_ldos.jld2.jld2")
xs, ys, zs = dict["xs"], dict["ys"], dict["mjdict"]
params, model = dict["params"], dict["model"]
plot_and_export(xs, ys, zs; filename =  "Output/FigEffective/FigEffective_test", labelbar = "dI/dV [G₀]", zfactor = 0.9, xlabel = "Φ/Φ₀", pages = false, layers = false)