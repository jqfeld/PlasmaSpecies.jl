using Revise
using PlasmaSpecies, Catalyst
using OrdinaryDiffEqDefault
using Symbolics
using QuadGK
using GLMakie
using AbstractTrees


# molecular constants for later use
energy(v::Int; ωe = 2_358.57, ωeχe = 14.324) = ωe * (v + 0.5) - ωeχe * (v + 0.5)^2 - ωe * 0.5 + ωeχe * 0.5^2
wavenumber_to_eV(wn) = wn / 8_065.54
Bv(v, Be = 1.99824, αe = 0.017318, γe = 0.0) = Be - αe * (v + 0.5) + γe * (v + 0.5)^2

const max_v = 45

F(Y) = Y == 0.0 ? 0.0 : Y^2 * quadgk(x -> exp(-x) / sinh(Y / sqrt(x))^2, 0, Inf)[1]
ΔE(v, w) = energy(v) + energy(w - 1) - energy(v - 1) - energy(w)
L = 0.2e-10
ħ = 1.0e-34
μ = 28 * 1.66e-27 / 2
kB = 1.38e-23
h = 2π * ħ
ω = 2_358.57 * 2 * pi * 3.0e10

TLT = π^2 * L^2 * ω^2 * μ / 2 / kB

Y(v, w, Tg) = sqrt(TLT / Tg) * abs(ΔE(v, w)) / 2_358.57

Y(2, 1, 300)

tr = SpeciesTree(
    [
        "N2[X,vib=$v]" for v in 0:max_v
    ]
)


##

p = 1.0e5
Tgas = 300.0
Ngas = p / kB / Tgas


energy_func(arg) = energy(parse(Int, arg.vib)) * h * 3.0e10

apply_energy!(tr[p"N2[X]"], energy_func)


##

rxs = PlasmaReaction[]

k0 = 4.0e-20

for v in 1:max_v, w in 1:max_v
    push!(
        rxs,
        PlasmaReaction(
            k0 * w * v * F(Y(v, w, Tgas)) * exp(ΔE(v, w) * h * 3.0e10 / 2 / kB / Tgas),
            ReactionFormula("N2[X,vib=$v] + N2[X,vib=$(w - 1)] --> N2[X,vib=$(v - 1)] + N2[X,vib=$(w)]")
        )
    )
end


reaction_energy.(rxs)

reaction_energy(rxs[1])

##
t = default_t()

@named sys = ReactionSystem(to_catalyst(tr, rxs...), t)

species(sys)

csys = complete(sys)

Nv1 = Ngas * 0.5 * kB * Tgas / (wavenumber_to_eV(energy(1)) * 1.6e-19)

prob = ODEProblem(
    csys,
    [
        csys.var"N2[X,vib=0]" => (Ngas - Nv1),
        csys.var"N2[X,vib=1]" => Nv1,
    ], (0.0, 1.0e-3)
)

sol = solve(prob)

##

sum(sol[:, end])

Ngas

##

scatter(
    [1;0;2:45], sol[:, end],
    axis = (yscale = log10,)
)

# lines!(
#     0 .. 10, v -> sol[2, end] * exp(-wavenumber_to_eV(energy(v)) * 1.6e-19 / kB / Tgas)
# )
