#region ## DEPENDENCES AND PARAMETERS
# using MKL

using Parameters, Arpack, ProgressMeter, LinearAlgebra, SparseArrays, Interpolations, Quantica
using Quantica: GreensFunction, allsitepositions, nsites, orbitals

@with_kw struct Params @deftype Float64  # nm, meV
    ħ2ome = 76.1996
    m0 = 0.023
    a0 = 10
    t = ħ2ome/(2m0*a0^2)
    Rcore = 60
    Rshell = Rcore
    Raxis = 0
    Vmax = 0
    Vmin = -40
    Vmax_lead = Vmax
    Vmin_lead = Vmin
    Rav = 0
    Δ0::ComplexF64 = 0.2
    g = 12
    μBΦ0 = 119.6941183 # meV nm^2
    P = 919.7 # meV nm
    Δg = 417  # meV
    Δs = 390  # meV
    α = P^2/3 * (1/Δg^2 - 1/(Δg + Δs)^2)  # nm^2, around 1.19 default
    α0 = 0
    Vexponent = 3
    echarge = 1
    L = 500 # the length, in nm, of the closed system
    δL = 0
    δR = 0
    δSL = δL/2
    δSR = δR/2
    δVL = 0
    δVR = 0
    ξd = 60
end

const σ0τx = @SMatrix[0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]
const σ0τy = @SMatrix[0 0 -im 0; 0 0 0 -im; im 0 0 0; 0 im 0 0]
const σ0τz = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]
const σ0τ0 = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
const σzτ0 = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 -1]
const σzτz = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
const σyτy = @SMatrix[0 0 0 -1; 0 0 1 0; 0 1 0 0; -1 0 0 0]
const σyτz = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 im; 0 0 -im 0]
const σyτ0 = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 -im; 0 0 im 0]
const σxτz = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 -1 0]
const σxτ0 = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]
const σ0 = SA[1 0; 0 1]
const σx = SA[0 1; 1 0]
const σy = SA[0 -im; im 0]
const σz = SA[1 0; 0 -1];

# Little-Parks and shell self-energy

function pairbreaking(n, Δ0, ξd, Rcore, Rshell)
    RLP = (Rcore + Rshell) / 2
    dshell = Rshell - Rcore
    Λ = ξd^2 * Δ0 / (1.76 * π * RLP^2) * (4 * (n - round(n))^2 + dshell^2 / RLP^2 * (n^2 + (round(n)^2)/3))
    #Δd = Δ0 * (1 - π/4 * Λ/Δ0 - π^2/32 * (Λ/Δ0)^2 - π^3/96 * (Λ/Δ0)^3)
    return Λ
end

function uUsadel(Δ0, Λ, ω)
    Δd = Δ0 * (1 - π/4 * Λ/Δ0 - π^2/32 * (Λ/Δ0)^2 - π^3/96 * (Λ/Δ0)^3)
    pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 
    6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
    nun = ω^2 - Δd^2 + Λ^2
    rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
    usa = 1/(2 * Δd) *(ω + sign(ω) * rai - sign(ω) * sqrt(complex(2 * ω^2 - 4 * nun/3 - nun^2 /(3 * pep) - pep/ 3 - sign(ω) * 2 * (Δd^2 + Λ^2) * ω / rai)))
    return usa
end


# function uUsadel(Δ0, Λ, ω)
#     Δd = Δ0 * (1 - π/4 * Λ/Δ0 - π^2/32 * (Λ/Δ0)^2 - π^3/96 * (Λ/Δ0)^3)
#     gap = complex(Δd^(2/3) - Λ^(2/3))^(3/2)
#     usa = ω / gap
#     return usa
# end

ΣS3DUsadel(Δ0, Λ, ω) = -(uUsadel(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1-uUsadel(Δ0, Λ, ω)^2))
#ΣS3DUsadel(Δ0, Λ, ω) = -im*(uUsadel(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(uUsadel(Δ0, Λ, ω)^2 - 1))
# ΣS3DUsadel(Δ0, Λ, ω) = -(uUsadel(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1- uUsadel(Δ0, Λ, ω)^2))

# ΣS3D(Δ, ω) = -(ω * σ0τ0 + Δ * σ0τx) / sqrt(Δ^2 - complex(ω)^2)

#endregion

#region ## CYLINDER MODEL, DECOMPOSED IN mⱼ

build_cyl(; nforced = nothing, kw...) = build_cyl(Params(; kw...); nforced)

function build_cyl(p::Params; nforced = nothing)
    @unpack Rcore, Rshell, Raxis, Rav, Vmin, Vmax, Vmin_lead, Vmax_lead, Vexponent,
            Δ0, ξd, g, μBΦ0, echarge, L, δL, δR, δSL, δSR, δVL, δVR, α, α0,
            ħ2ome, t, a0, m0 = p

    R´ = floor(Rcore/a0)*a0
    L´ = L + δSL + δSR

    lat = if iszero(Rav)
            LP.square(; a0, names = :C) |> unitcell((1, 0), region = r -> max(a0, Raxis) <= r[2] <= R´) |>
            transform!(r -> r + SA[(δL + L´- δR)/2, 0]) # center of island
        else
            LP.square(; a0, names = :C) |> unitcell((1, 0)) |> transform!(r -> r + SA[(δL + L´- δR)/2, Rav]) # center of island
        end

    # The radial p² is derived in https://doi.org/10.1016/j.cpc.2015.08.002 (Eq. (20) with l = 0 and constant m)
    # The resulting eigenstates are F(r) = Ψ(r) * √r
    # For constant mass it reads t * (2Fⱼ - Fⱼ₊₁ * ½(rⱼ + rⱼ₊₁)/√(rⱼrⱼ₊₁) - Fⱼ₋₁ * ½(rⱼ + rⱼ₋₁)/√(rⱼrⱼ₋₁)) =
    # = onsite(2t) - hopping((r, dr) -> t*r/√(r+dr/2)(r-dr/2)) where t = ħ²/(2ma₀²) and Fⱼ=Ψ(rⱼ)√rⱼ
    # Boundary condition: remove j = 0 and replace the j = 1 onsite by 3t/2
    #p² = onsite(r -> σ0τz * of(r)) + hopping((r, dr) -> tf(r, dr) * σ0τz, range = a0)
    p² = onsite(r -> σ0τz * t * ifelse(r[2] ≈ a0, 2.0 + 1.5, 2.0 + 2.0)) -
        hopping((r, dr) -> t * σ0τz * ifelse(iszero(dr[1]), r[2]/sqrt(r[2]^2 - 0.25*dr[2]^2), 1), range = a0)

    V(ρ, v0, v1) =  v0 + (v1-v0) * (ρ/R´)^Vexponent
    dϕ(ρ, v0, v1) = -(Vexponent / R´) * (v1-v0) * (ρ/R´)^(Vexponent-1) # ϕ=-V
    V(r) = smooth(r[1], (0, δSL, L´-δSR, L´), (V(r[2], Vmax_lead, Vmin_lead), V(r[2], Vmax, Vmin), V(r[2], Vmax, Vmin))) +
           barrier(r[1], (0, δL), δVL) + barrier(r[1], (L´-δR, L´), δVR)
    dϕ(r) = smooth(r[1], (0, δL, L´-δR, L´), (dϕ(r[2], Vmax_lead, Vmin_lead), dϕ(r[2], Vmax, Vmin), dϕ(r[2], Vmax, Vmin)))
    potential = onsite(r -> V(r) * σ0τz)
    # <z|pz|psi> = <z|-i∂z|psi>  = -i[psi(z+a0)-psi(z-a0)]/(2a0) = i dr[1]/(2a0^2)
    rashba = hopping((r, dr) -> (α0 + α * dϕ(r)) * (im * dr[1] / (2a0^2)) * σyτz)

    area_LP = pi * (0.5*Rcore + 0.5*Rshell)^2
    # eA/ħ = 1/2 * pi * n r/area_LP, where n = area_LP*B/Φ0 and Φ0 = h/2e
    # V_Zeeman = 1/2 g μB Φ0 n/area_LP
    zeeman! = @onsite!((o; n = 0) -> o + (0.5 * g * μBΦ0 * n / area_LP) * σzτ0)
    eAφ(r, n) = echarge * 0.5 * pi * n * r[2]/ area_LP
    nint(n) = ifelse(nforced === nothing, round(Int, n), nforced)
    mj(Z, n) = Z + ifelse(iseven(nint(n)), 0.5, 0.0)
    J(Z, n) = mj(Z, n)*σ0τ0 - 0.5*σzτ0 - 0.5*nint(n)*σ0τz
    gauge! = @onsite!((o, r; n = 0, Z = 0) -> o +
        (eAφ(r, n) * σ0τz + J(Z, n) / r[2])^2 * σ0τz * ħ2ome / (2 * m0) -
        (eAφ(r, n) * σ0τz + J(Z, n) / r[2]) * (α0 + α * dϕ(r[2], Vmax, Vmin)) * σzτz)

    # Little-Parks and self-energy from the shell
    Λ(n) = pairbreaking(n, Δ0, ξd, Rcore, Rshell)
    τprofile(r) = smooth(r[1], (δSL, δSL, L´-δSR, L´-δSR), (0, 1.0, 1.0))
    ΣS! = @onsite!((o, r; ω = 0, n = 0, τΓ = 0.0) ->
        o + τprofile(r) * (τΓ * Δ0 * ΣS3DUsadel(Δ0, Λ(n), ω));
        # o + τprofile(r) * (τΓ * Δ0 * ΣS3D(Δ0, ω));
        region = iszero(Rav) ? r -> r[2] > R´ - a0/2 : Returns(true))

    #OldΣS! = @onsite!((o, r; ω = 0, τΓ = 1, Γ = Δ, n = 0, Δ0 = 0.0, slead = false) ->
    #    o - Δ0 * littleparks(n) * σ0τx + Δsmooth(r, slead) * ΣS(τΓ, Γ, littleparks(n) * Δ, ω); region = r -> r[2] >= R´ - 0.1 * a0)

    # Leads
    latlead = lat

    potential_lead = onsite(r -> V(r[2], Vmax_lead, Vmin_lead) * σ0τz)
    hlead = hamiltonian(latlead, p² + potential_lead + rashba, orbitals = Val(4))
    # phlead = slead ? hlead |> parametric(gauge!, zeeman!, ΣS_lead!, ΣΔcore_lead!) : hlead |> parametric(gauge!, zeeman!)
    phNlead = hlead |> parametric(gauge!, zeeman!)
    phSlead = hlead |> parametric(gauge!, zeeman!, ΣS!)

    # same as above, but use Vmax/min instead of Vmax_lead/min_lead
    potential_lead´ = onsite(r -> V(r[2], Vmax, Vmin) * σ0τz)
    hlead´ = hamiltonian(latlead, p² + potential_lead´ + rashba, orbitals = Val(4))
    phSlead´ = hlead´ |> parametric(gauge!, zeeman!, ΣS!)

    # Normal Lead self-energy
    function ΣNlead(ω, side = 1; kw...)
        h´ = phNlead(; ω, kw...)
        g0 = greens(h´, Schur1D(), boundaries = (0,))(ω, side=>side)
        return h´[(-side,)] * g0 * h´[(side,)]
    end

    # S lead self-energy
    function ΣSlead(ω, side = 1; kw...)
        h´ = phSlead(; ω, kw...)
        g0 = greens(h´, Schur1D(), boundaries = (0,))(ω, side=>side)
        return h´[(-side,)] * g0 * h´[(side,)]
    end

    img(Σ) = (Σ'-Σ)/2im

    # Central
    # latL = LP.square(; a0) |> unitcell(region = r -> 0 <= r[1] <= L´ && Raxis < r[2] <= R´)
    latL = latlead |> unitcell(region = r -> 0 <= r[1] <= L´)
    iC = siteindices(latL, region = r -> L´/2 <= r[1] <= L´/2+0.99a0) |> collect
    # iLL = siteindices(latL, region = r -> δL <= r[1] <= δL+0.99a0) |> collect
    iL = siteindices(latL, region = r -> 0 <= r[1] < 0.99a0) |> collect
    iR = siteindices(latL, region = r -> L´-0.99a0 < r[1] <= L´) |> collect
    # iS = siteindices(latL, region = r -> r[2] > R´ - 0.99*a0) |> collect
    ΣN0 = onsite(0I, indices = vcat(iL, iR)) + hopping(0I, range = Inf, indices = (iR => iR, iL => iL))
    h0 = hamiltonian(latL, p² + potential + rashba + ΣN0, orbitals = Val(4))

    ΣR! = @block!((b; τR = 1.0, ω = 0, kw...) -> iszero(τR) ? b : b + τR * ΣSlead(ω, 1; kw...), iR)
    ΣL! = @block!((b; τL = 1.0, ω = 0, kw...) -> iszero(τL) ? b : b + τL * ΣNlead(ω, -1; kw...), iL)

    ph = h0 |> parametric(gauge!, zeeman!, ΣS!, ΣR!, ΣL!)

    hΓ = hamiltonian(latL, ΣN0, orbitals = Val(4))
    ΓL! = @block!((b; τL = 1, ω = 0, kw...) -> iszero(τL) ? b : b + τL * img(ΣNlead(ω, -1; kw...)), iL)
    ΓL = (hΓ |> parametric(ΓL!))[:, iL]


    return ph, iC, iL, ΓL, phSlead´
end

function smooth(x, (x0, x1)::NTuple{2,Any}, (y0, y1)::NTuple{2,Any})
    x <= x0 ? y0 :
    x >= x1 ? y1 :
    y0 + (y1-y0) * (tanh(5*(x-0.5*(x0+x1))/abs(x1-x0))+1)/2
end

smooth(x, (x0, x1, x2)::NTuple{3,Any}, (y0, y1, y2)::NTuple{3,Any}) =
    smooth(x, (x0, x1), (y0, y1)) + smooth(x, (x1, x2), (0, y2-y1))
smooth(x, (x0, x1, x2, x3)::NTuple{4,Any}, (y0, y1, y2)::NTuple{3,Any}) =
    smooth(x, (x0, x1), (y0, y1)) + smooth(x, (x2, x3), (0, y2-y1))

barrier(x, (x0, x1), V) = x0 < x < x1 ? V*exp(-25*(x-0.5*(x0+x1))^2/(x1-x0)^2) : 0.0

#endregion

#region ## HEXAGONAL MODEL IN 3D - SYMMETRY-REDUCED

function tohex((x, y, z))
    z = x + y*im
    φ = mod(angle(z), pi/3)
    return sqrt(x*x + y*y) * (cos(φ) + sin(φ)/√3)
end

function dxtohex((x, y, z))
    abscosϕ2 = x*x / (x*x + y*y)
    dx = sign(x) * ifelse(abscosϕ2 < 0.25, 0, 1)
    return dx
end

function dytohex((x, y, z))
    abssinϕ2 = y*y / (x*x + y*y)
    dy = sign(y) * ifelse(abssinϕ2 < 0.75, 1, 2)/√3
    return dy
end

build_hex(; kw...) = build_hex(Params(; kw...))

function build_hex(p::Params)
    @unpack Rcore, Rshell, Raxis, Vmin, Vmax, Vmin_lead, Vmax_lead, Vexponent,
            Δ0, γ, maxlobe2, g, μBΦ0, echarge, L, δL, δR, δSL, δSR, δVL, δVR, α, α0,
            ħ2ome, a0, m0, as, ms, Vs, Ls = p

    δh = 1 - γ/4   # normalized gap at \Phi = 0.5 * \Phi_0
    δ2 = maxlobe2

    R´ = a0/6 + floor((Rcore-a0/2)/a0)*a0 + a0/2

    ts = ħ2ome/(2*ms*as^2)
    t = 2*ħ2ome/(m0*a0^2)
    acc = a0/√3
    t1 = t/sqrt((1+2as)/(3acc))
    t2 = ts/sqrt(0.5+0.25*acc/as)
    o = 3t
    o1 = 2*ħ2ome/(2as+acc) * (1/(ms*as) + 2/(m0*acc))
    tz = ħ2ome/(2m0*a0^2)
    oz = 2tz

    ϕr(r) = atan(r[2], r[1])
    dϕr(r, dr) = ϕr(r+dr/2) -  ϕr(r-dr/2)

    lat0 = LP.honeycomb(; a0, dim = Val(3), bravais = (a0 .* (1/2, √3/2, 0), a0 .* (-1/2, √3/2, 0), a0 .* (0, 0, 1)))
    transform!(r -> r + SA[a0/2, 0, 0], lat0)
    irred(r) = 0 <= ϕr(r) < 0.9999 * pi/3
    latC = lat0 |> unitcell((0, 0, 1), region = r -> tohex(r/Raxis) > 1 && tohex(r/R´) < 0.9999 && irred(r))|> lattice(names = (:A,:B))
    latI = lat0 |> unitcell((0, 0, 1), region = r -> 0.9999 <= tohex(r/R´) <= 1.0001 && irred(r)) |> lattice(names = (:AS,:BS))
    lat = combine(latC, latI)

    p² = onsite((o + oz)*σ0τz, sublats = (:A, :B)) +
         onsite((o1 + oz)*σ0τz, sublats = (:AS, :BS)) +
         hopping(-t*σ0τz, sublats = (:A,:B) .=> (:B,:A), range = acc) +
         hopping(-t1*σ0τz, sublats = (:A,:B,:AS,:BS,:AS,:BS) .=> (:BS,:AS,:B,:A,:BS,:AS), range = acc) +
         hopping((r, dr) -> iszero(dr[3]) ? 0*σ0τz : -tz*σ0τz, range = a0)

    rot = [cos(pi/3) -sin(pi/3) 0; sin(pi/3) cos(pi/3) 0; 0 0 1]
    rdr(r1, r2) = (0.5*(r1+r2), r2-r1)
    rdrwrap(r, dr) = argmin(((r, dr),) -> norm(dr), ((r, dr), rdr((r - dr/2), rot * (r + dr/2)), rdr((r - dr/2), rot' * (r + dr/2))))
    reg(r, dr) = norm(last(rdrwrap(r, dr))) <= 1.001 * acc
    p² += hopping(-t*σ0τz, sublats = (:A,:B) .=> (:A,:B), range = 2R´, region = reg) +
          hopping(-t1*σ0τz, sublats = (:A,:B,:AS,:BS,:AS,:BS) .=> (:AS,:BS,:A,:B,:AS,:BS), range = 2R´, region = reg)

    potential(r) = (Vmax + (Vmin-Vmax) * tohex(r/R´)^Vexponent) * σ0τz
    rashba(r, dr) = α/2a0 * Vexponent * (Vmin-Vmax)/R´ * tohex(r/R´)^(Vexponent-1) *
        (-im * dxtohex(r/R´) * (dr[3]*σyτz - dr[2]*σzτ0) - im * dytohex(r/R´) * (dr[1]*σzτ0 - dr[3]*σxτ0))

    h = lat |> hamiltonian(p² + onsite(potential) + hopping(rashba), orbitals = Val(4))

    # eA/ħ = 1/2 * pi * n r/area_LP, where n = area_LP*B/Φ0 and Φ0 = h/2e
    # V_Zeeman = 1/2 g μB Φ0 n/area_LP
    area_LP = pi * (0.5*Rcore + 0.5*Rshell)^2
    eA(r, n) = echarge * 0.5 * pi * n * SA[r[2], -r[1], 0]/ area_LP
    nambudiag(t) = Diagonal(SA[t, t, conj(t), conj(t)])
    Zn(Z, n) = Z + 0.5 * mod(round(n), 2)  ## For SU(2) even/odd boundary conditions

    peierls(n, (r, dr)) = dr' * eA(r, n) + 0.5 * round(n) * dϕr(r, dr)

    flux! = @hopping!((t, r, dr; n = 0) -> t * nambudiag(cis(peierls(n, rdrwrap(r, dr)))))
    flux0! = @onsite!((o, r; n = 0, Z = 0) ->
        o - t * σ0τz * 2 * real.(nambudiag(cis(peierls(n, rdr(r, rot*r)))) * cispi(Zn(Z, n)/3));
        region = r -> norm(r) < 1.001*acc)
    zeeman! = @onsite!((o; n = 0) -> o + (0.5 * g * μBΦ0 * n / area_LP) * σzτ0)
    wrap! = @hopping!((t, r, dr; n = 0, Z = 0) ->
        norm(dr) > 1.001*acc ? t*cispi(-sign(dϕr(r, dr))*(Zn(Z, n))/3) : t)

    # Shell self-energy
    # Δmin = 0.00001
    # littleparks(n) = max(Δmin, 1 - γ * (n - round(n))^2 - 0.25*(1 - maxlobe2)*n^2)
    lp(n) = littleparks(n, δh, δ2)
    phShell = LP.linear(; a0 = as) |> hamiltonian(onsite((2ts + Vs) * σ0τz) - hopping(ts*σ0τz), orbitals = Val(4)) |>
              parametric(@onsite!((o; Δshell = Δ0) -> o + Δshell * σ0τx))

    function ΣS(ω, Δ0)
        h´ = phShell(; Δshell = Δ0)
        if isinf(Ls)
            g0 = greens(h´, Schur1D(), boundaries = (0,))(ω, 1=>1)
            Σ = (σ0τz * g0[1,1] * σ0τz) * t2^2
        else
            L´ = round(Int, Ls/as)
            os = orbitalstructure(h´)
            g₀ω = greens(h´, Schur1D(), boundaries = (0,))(ω)
            g₀_L´L´⁻¹ = Quantica.unflatten_blocks(inv(flatten(g₀ω[L´=>L´], os)), os)
            g11 = g₀ω[1=>1] - g₀ω[L´=>1] * g₀_L´L´⁻¹ * g₀ω[1=>L´]
            Σ = (σ0τz * g11[1,1] * σ0τz) * t2^2
        end
        return Σ
    end

    ΣS! = @onsite!((o; ω = 0, n = 0, τ = 1.0) ->
        o + τ * ΣS(ω, Δ0 * lp(n)); sublats = (:AS, :BS))

    ph = h |> parametric(flux!, flux0!, wrap!, zeeman!, ΣS!)

    return ph
end

#endregion

#region ## DOS

## Finite length, using greens
function dos(hf::Function, xrng, ωrng, L, ph)
    os = orbitalstructure(parent(ph))
    pts = Iterators.product(xrng, ωrng)
    d = @showprogress pmap(pts) do pt
        x, ω = pt
        h = hf(x, ω)
        g₀ = greens(h, Schur1D(), boundaries = (0,))
        g₀ω = g₀(ω)
        L´ = round(Int,L + 1)
        g₀_L´L´⁻¹ = Quantica.unflatten_blocks(inv(flatten(g₀ω[L´=>L´], os)), os)
        g11 = g₀ω[1=>1] - g₀ω[L´=>1] * g₀_L´L´⁻¹ * g₀ω[1=>L´]
        return -imag(tr(tr(g11)))/π
    end
    return reshape(d, size(pts)...)
end

## Semi-infinite, using greens
function dos(hf::Function, xrng, ωrng; kw...)
    pts = Iterators.product(xrng, ωrng)
    d = @showprogress pmap(pts) do pt
        x, ω = pt
        gf = greens(hf(x, ω), Schur1D(), boundaries = (0,))
        return -imag(tr(tr(gf(ω, 1=>1))))/π
    end
    return reshape(d, size(pts)...)
end

function dos_batch(hf::Function, xrng, ωrng; kw...)
    pts = Iterators.product(xrng, ωrng)
    d = @showprogress @distributed vcat for pt in collect(pts)
        x, ω = pt
        gf = greens(hf(x, ω), Schur1D(), boundaries = (0,))
        dos = -imag(tr(tr(gf(ω, 1=>1))))/π
        dos
    end
    return reshape(fetch(d), size(pts)...)
end

## Semi-infinite, using greens
function dos(hf::Function, xrng, ωrng; kw...)
    pts = Iterators.product(xrng, ωrng)
    d = @showprogress pmap(pts) do pt
        x, ω = pt
        gf = greens(hf(x, ω), Schur1D(), boundaries = (0,))
        return -imag(tr(tr(gf(ω, 1=>1))))/π
    end
    return reshape(d, size(pts)...)
end

function dos_r(hf::Function, x, ωrng; kw...)
    pts = ωrng
    d = @showprogress pmap(pts) do pt
        ω = pt
        gf = greens(hf(x, ω), Schur1D(), boundaries = (0,))
        return -imag.(tr.(diag(gf(ω, 1=>1)))) ./ π
    end
    return hcat(d...)
end

## Finite length, using \ (ldiv)
function dos_ldiv(hf::Function, xrng, ωrng, iL, ph)
    IL = identity_inds(iL, size(ph, 1))
    pts = Iterators.product(xrng, ωrng)
    x, ω = first(pts)
    m = similarmatrix(hf(x, ω), flatten)
    d = @showprogress pmap(pts) do pt
        x, ω = pt
        h = hf(x, ω)
        bloch!(m, h)
        IL´ = (ω*I - m) \ IL
        dos = 0.0
        for (c, i) in enumerate(iL), k in 1:4
            dos += -imag(IL´[4(i-1)+k, 4(c-1)+k])/π
        end
        return dos
    end
    return reshape(d, size(pts)...)
end

function identity_inds(iL, n)
    IL = fill(0.0, 4*n, 4*length(iL))
    for (c, i) in enumerate(iL), k in 1:4
        IL[4(i-1) + k, 4(c-1) + k] = 1
    end
    return IL
end

#endregion

#region ## CONDUCTANCE

# function conductance(ph, Γ, nrng, ωrng; kw...)
#     pts = Iterators.product(nrng, ωrng)
#     hmat = similarmatrix(ph, flatten)
#     Γmat = similarmatrix(Γ, Matrix{Quantica.blockeltype(Γ)})
#     Γmat´ = similar(Γmat)
#     iΓ = Quantica.axesflat(Γ, 2)
#     τe, τz = SA[1,1,0,0],  SA[1,1,-1,-1]
#     gs = @showprogress pmap(pts) do pt
#         n, ω = pt
#         cond = 0.0
#         bloch!(hmat, ph(; kw..., ω=ω, n=n))
#         bloch!(Γmat, Γ(; kw..., ω=ω, n=n))
#         copy!(Γmat´, Γmat)

#         # Why is this faster than a straight lu(hmat) ???
#         luig = lu(ω*I - hmat)
#         GrΓ = view(ldiv!(luig, Γmat), iΓ, :)
#         GaΓ = view(ldiv!(luig', Γmat´), iΓ, :)

#         # G = 2i*Tr[(GrΓ - GaΓ)τe] + 4*Tr[GaΓ τh GrΓ τe] - 4*Tr[GaΓ τe GrΓ τe]
#         #   = 2i*Tr[(GrΓ - GaΓ)τe] - 4*Tr[GaΓ τz GrΓ τe]
#         cond -= 4 * real(trace_product(GaΓ, τz, GrΓ, τe))
#         GrΓ .-= GaΓ
#         cond += 2 * real(im * trace_product(GrΓ, τe))

#         return cond # The unit is e^2/h
#     end
#     return gs
# end

function conductance(ph, Γ, nrng, ωrng; kw...)
    pts = Iterators.product(nrng, ωrng)
    hmat = similarmatrix(ph, flatten)
    Γmat = similarmatrix(Γ, Matrix{Quantica.blockeltype(Γ)})
    Γmat´ = similar(Γmat)
    iΓ = Quantica.axesflat(Γ, 2)
    τe, τz = SA[1,1,0,0],  SA[1,1,-1,-1]

    # This is better than pmap because pmap sends all closed-over arrays at each iteration,
    # while @distributed does only once per batch
    gs = @showprogress @distributed vcat for pt in collect(pts)
        n, ω = pt
        cond = 0.0
        bloch!(hmat, ph(; kw..., ω=ω, n=n))
        bloch!(Γmat, Γ(; kw..., ω=ω, n=n))
        copy!(Γmat´, Γmat)

        # Why is this faster than a straight lu(hmat) ???
        luig = lu(ω*I - hmat)
        GrΓ = view(ldiv!(luig, Γmat), iΓ, :)
        GaΓ = view(ldiv!(luig', Γmat´), iΓ, :)

        # G = 2i*Tr[(GrΓ - GaΓ)τe] + 4*Tr[GaΓ τh GrΓ τe] - 4*Tr[GaΓ τe GrΓ τe]
        #   = 2i*Tr[(GrΓ - GaΓ)τe] - 4*Tr[GaΓ τz GrΓ τe]
        cond -= 4 * real(trace_product(GaΓ, τz, GrΓ, τe))
        GrΓ .-= GaΓ
        cond += 2 * real(im * trace_product(GrΓ, τe))
        cond  # The unit is e^2/h
    end
    return reshape(fetch(gs), size(pts)...)
end

function trace_product(GΓ::AbstractMatrix{T}, τ::SVector{N}) where {N,T}
    tp = zero(real(T))
    cols = size(GΓ, 2)
    for j in 1:cols
        tp += GΓ[j,j] * τ[mod1(j, N)]
    end
    return tp
end

function trace_product(GΓ1::AbstractMatrix{T}, τ1::SVector{N}, GΓ2, τ2::SVector{N}) where {N,T}
    rows, cols = size(GΓ1)
    tr = zero(real(T))
    for j in 1:cols, i in 1:rows
        tr += GΓ1[i, j] * τ1[mod1(j, N)] * GΓ2[j, i] * τ2[mod1(i, N)]
    end
    return tr
end

#endregion

#region ## PLOTTING AND EXPORTING

using VegaLite, CairoMakie, JLD2

function plotdata(dos::AbstractMatrix;
        xmax = 3, xlims = (0, xmax), ymax = 0.22, ylims = (0, ymax), zmin = 0, zmax = Inf, zlims = (zmin, zmax),
        xlabel = "n", ylabel = "V [mV]", labelbar = "DOS [a.u.]", aspect = 1.6, label = "", kw...)
    xs = range(xlims..., length = size(dos, 1))
    ys = range(ylims..., length = size(dos, 2))
    fig = Figure(resolution = (1200, round(Int, 1200/aspect)), font =:sans)
    ax = Axis(fig; aspect, xlabel, ylabel)
    dosc = clamp.(dos, zmin, zmax)
    hmap = CairoMakie.heatmap!(xs, ys, dosc; colormap = :thermal, colorrange = zlims, kw...)
    cbar = Colorbar(fig, hmap, label = labelbar,
        labelpadding = 5,
        ticklabelpad = 30)
    fig[1, 1] = ax
    fig[1, 2] = cbar
    Label(fig[0, :], string(label), fontsize = 20)
    return fig
end

function plot_and_export(xs, ys, mjdict; filename, zfactor = 0.3, zfactor´ = zfactor, xlabel = "Γ/t", layers = true, kw...)
    path = dirname(filename)
    path0 = pwd()

    try run(`mkdir -p $path`) catch end
    @save(filename*".jld2", xs, ys, mjdict, params, model, filename)

    data = values(mjdict)
    ms = keys(mjdict)
    p = plotdata(sum(data); xlims = extrema(xs), zmax = zfactor * maximum(sum(data)), 
                 ylims = extrema(real.(ys)), xlabel, label = "total", kw...);

    if !layers
        save("$filename.pdf", p)
    else
        filelist = "0.tmp.pdf"
        save("$path/0.tmp.pdf", p);

        counter = 0
        for (mj, d) in sort!(collect(zip(ms, data)))
            pn = plotdata(d; xlims = extrema(xs), zmax = zfactor´ * maximum(d),
                          ylims = extrema(real.(ys)), xlabel, label = "mj = $mj", kw...);
            counter += 1; save("$path/$counter.tmp.pdf", pn); filelist *= " $counter.tmp.pdf"
        end
        cd(path)
        try
            pdftk = "pdftk $filelist cat output out.pdf"
            run(`bash -c $pdftk`)
            run(`mv out.pdf $(basename(filename)).pdf`)
        catch
            throw(error())
        finally
            c = "rm *.tmp.pdf"
            run(`bash -c $c`)
            cd(path0)
        end
    end
    return p
end

#endregion