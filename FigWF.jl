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

#region ## FUNCTIONS

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

function figWF_Raxis(path; n = 0, E0 = 0, τΓ = 0, phi = 0)
    try run(`mkdir -p $path`) catch end
    f = Figure(resolution = (300, 600), font = "CMU Serif Roman")
    Axis(f[1, 1], xticks = 0:10:70, yticks = -0.4:0.05:0.4, xlabel = L"r(nm)", yticksvisible = false, yticklabelsvisible = false)
    cols = (:red, :orange, :green, :blue)
    for (i, Raxis) in enumerate(0:10:60)
        a0 = 1
        @show Raxis, a0
        Lsys, iC, iL, ΓL, leadS = build_cyl(; model..., a0, Raxis = Raxis + a0, Rcore = 70, Vmin = 0, Vmax = 0);
        for Z in 0:3
            hh = leadS(; Z, ω = 0, τΓ, n)
            ε, F = spectrum(hh, phi; method = ArpackPackage(; nev = 4, sigma = E0 + 0.1im))[around = E0]
            psi =  map(norm, F ./ sqrt.(last.(sitepositions(hh)))) |> vec
            @show ε
            iszero(Z) && iszero(Raxis) ? pushfirst!(psi, first(psi)) : pushfirst!(psi, 0)
            iszero(τΓ) && push!(psi, 0)
            psi .*= 1.2
            psi .+= 0.05 * (i-1)
            sp = last.(sitepositions(leadS()))
            rs = prepend!(sp, Raxis)
            iszero(τΓ) && push!(sp, 70)
            lines!(rs, psi, color = cols[Z+1], ylims = (-0.03, 0.43))
            band!(rs, 0.05 * (i-1), psi, color = (cols[Z+1], 0.1), ylims = (-0.03, 0.43))
        end
    end
    save("$(path).pdf", f)
    return f
end

function figWF_dome(path; n = 0, E0 = 0, τΓ = 0, phi = 0)
    try run(`mkdir -p $path`) catch end
    f = Figure(resolution = (300, 600), font = "CMU Serif Roman")
    Axis(f[1, 1], xticks = 0:10:70, yticks = -0.4:0.05:0.4, xlabel = L"r(nm)", yticksvisible = false, yticklabelsvisible = false)
    cols = (:red, :orange, :green, :blue)
    for (i, V) in enumerate(0:20:120)
        a0 = 1
        @show V
        Lsys, iC, iL, ΓL, leadS = build_cyl(; model..., a0, Rcore = 70, Vmin = 0, Vmax = V);
        for Z in 0:3
            hh = leadS(; Z, ω = 0, τΓ, n)
            ε, F = spectrum(hh, phi; method = ArpackPackage(; nev = 4, sigma = E0 + 0.1im))[around = E0]
            psi =  map(norm, F ./ sqrt.(last.(sitepositions(hh)))) |> vec
            @show ε
            iszero(Z) ? pushfirst!(psi, first(psi)) : pushfirst!(psi, 0)
            iszero(τΓ) && push!(psi, 0)
            psi .*= 1.2
            psi .+= 0.05 * (i-1)
            sp = last.(sitepositions(leadS()))
            rs = prepend!(sp, 0)
            iszero(τΓ) && push!(sp, 70)
            lines!(rs, psi, color = cols[Z+1], ylims = (-0.03, 0.43))
            band!(rs, 0.05 * (i-1), psi, color = (cols[Z+1], 0.1), ylims = (-0.03, 0.43))
        end
    end
    save("$(path).pdf", f)
    return f
end


#endregion

#region ## CONFIGURATION FIG 5 - COMMON

base_model = (;
    a0 = 2,
    ξd = 70, Rshell = 70, Rcore = 70, Raxis = 50,
    Δ0 = 0.23, α=0, g = 0, Vexponent = 2)

base_params = (; τΓ = 0.6)

ωrng = range(-0.26, 0.26, 300) .+ 0.001im
nrng = range(0, 2.499, 400)

# ωrng = range(-0.26, 0.26, 60) .+ 0.001im
# nrng = range(2, 2.49, 80)

function density_of_states()
    xs, ys = nrng, ωrng;
    mjdict = Dict([(@show Z; Z => dos(nω_leadS(Z), xs, ys)) for Z in mjs]);
    return xs, ys, mjdict
end

function plot((xs, ys, mjdict), filename;  zfactor = 1.0, zfactor´ = 1.0, xlabel = "Φ/Φ₀", labelbar = "LDOS [a.u]")
    plot_and_export(xs, ys, mjdict; filename, zfactor, zfactor´, xlabel, labelbar)
end

#endregion


#region ## PANEL 5(a)

model = (; base_model..., a0 = 2, Raxis = 50, Vmin = -34, Vmax = -34, ξd = 0.001)
params = (; base_params..., τΓ = 250)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

data = density_of_states()
plot(data, "Output/FigWF/FigWFa"; zfactor = 0.9, zfactor´ = 0.6)

#endregion

#region ## PANEL 5(b)

model = (; base_model..., a0 = 2, Raxis = 50, Rshell = 90, Vmin = -34, Vmax = -34, ξd = 0.0001)
params = (; base_params..., τΓ = 250)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

data = density_of_states()
plot(data, "Output/FigWF/FigWFb"; zfactor = 0.9, zfactor´ = 0.6)

#endregion

#region ## PANEL 5(b) Inset

model = (; base_model..., a0 = 2, Raxis = 50, Rshell = 90, Vmin = -34, Vmax = -34, ξd = 0.0001)
params = (; base_params..., τΓ = 250)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; nforced = 1, model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

data = density_of_states()
plot(data, "Output/FigWF/FigWFbinset"; zfactor = 0.9, zfactor´ = 0.6)

#endregion

#region ## PANEL 5(c)
# ξ=140nm, Rshell = 70nm, dshell = 0
model = (; base_model..., a0 = 2, Raxis = 50, Vmin = -34, Vmax = -34, ξd = 140)
params = (; base_params..., τΓ = 250)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

data = density_of_states()
plot(data, "Output/FigWF/FigWFc"; zfactor = 0.9, zfactor´ = 0.6)

#endregion

#region ## PANEL 5(d)
# ξ=140nm, Rshell = 90nm, dshell = 20nm
model = (; base_model..., a0 = 2, Raxis = 50, Rshell = 90, Vmin = -34, Vmax = -34, ξd = 140)
params = (; base_params..., τΓ = 250)
mjs = -7:7

Lsys, iC, iL, ΓL, leadS = build_cyl(; model...);

nω_leadS(Z) = (n, ω)  -> leadS(; Z = Z, ω = ω, n, τΓ = params.τΓ)

data = density_of_states()
plot(data, "Output/FigWF/FigWFd"; zfactor = 0.5, zfactor´ = 0.6)

#endregion


#region ## FIG WF Raxis (e)

model = (;
    a0 = 1,
    ξd = 0.001, Rshell = 70, Rcore = 70, Raxis = 0,
    Δ0 = 0.23, α = 0, g = 0, Vexponent = 2)

base_params = (; τΓ = 250)

figWF_Raxis("Output/figWF/figWF_Raxis"; τΓ=0)

#endregion

#region ## FIG WF Dome (f)

model = (;
    a0 = 1,
    ξd = 0.001, Rshell = 70, Rcore = 70, Raxis = 0,
    Δ0 = 0.23, α = 0, g = 0, Vexponent = 2)

base_params = (; τΓ = 250)

figWF_dome("Output/figWF/figWF_dome_1"; n = 0, τΓ = 0)

#endregion



## Lattice
hh = leadS(; Z = 8, ω = 0+0.0001im, n = 0, τΓ = 1)
vlplot(hh, maxthickness = 0.06)

## Lattice
hh = Lsys(; Z = 8, ω = 0.0+0.0001im, n = 1, τΓ = 0.0)
vlplot(hh, maxthickness = 0.01)

## Bands
Lsys, iC, iL, ΓL, leadS = build_cyl(; model..., a0 = 5, Raxis = 0, Rcore = 70, Vmin = -120, Vmax = 0);

bands(0; nev = 10, τΓ = 1, n = 0)

## Load and plot
dict = load("Output/FigEffective/FigEffective_CZ_cond_long.jld2")
xs, ys, zs = dict["xs"], dict["ys"], dict["mjdict"]
params, model = dict["params"], dict["model"]
plot_and_export(xs, ys, zs; filename =  "Output/FigEffective/FigEffective_CZ_cond_long", labelbar = "dI/dV [G₀]", zfactor = 0.2, xlabel = "Φ/Φ₀", pages = false, layers = false)