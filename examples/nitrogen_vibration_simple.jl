using Revise
using PlasmaSpecies, Catalyst
using OrdinaryDiffEqDefault
using Symbolics
using QuadGK
using GLMakie
using AbstractTrees

t = default_t()
D = default_time_deriv()

# molecular constants for later use
energy(v::Int; ωe = 2_358.57, ωeχe = 14.324) = ωe * (v + 0.5) - ωeχe * (v + 0.5)^2 - ωe * 0.5 + ωeχe * 0.5^2
wavenumber_to_eV(wn) = wn / 8_065.54
Bv(v, Be = 1.99824, αe = 0.017318, γe = 0.0) = Be - αe * (v + 0.5) + γe * (v + 0.5)^2

const max_v = 1

ΔE(v) = energy(v) - energy(v - 1)

L = 0.2e-10
ħ = 1.0e-34
μ = 28 * 1.66e-27 / 2
kB = 1.38e-23
h = 2π * ħ
ω = 2_358.57 * 2 * pi * 3.0e10

tr = SpeciesTree(
    [
        "N2[X,vib=$v]" for v in 0:max_v
    ]
)


##


energy_func(arg) = energy(parse(Int, arg.vib)) * h * 3.0e10
apply_energy!(tr[p"N2[X]"], energy_func)

##

@variables Tgas(t)

rxs = PlasmaReaction[]

kVV0 = 1.0e-20
kVT0 = 3.0e-10

for v in 1:max_v
    for w in 0:max_v
        push!(
            rxs,
            PlasmaReaction(
                kVT0 * exp(ΔE(v) * h * 3.0e10 / 2 / kB / Tgas),
                ReactionFormula("N2[X,vib=$v] + N2[X,vib=$w] --> N2[X,vib=$(v - 1)] + N2[X,vib=$w]")
            )
        )
        push!(
            rxs,
            PlasmaReaction(
                kVT0 * exp(-ΔE(v) * h * 3.0e10 / 2 / kB / Tgas),
                ReactionFormula("N2[X,vib=$(v - 1)] + N2[X,vib=$w] --> N2[X,vib=$(v)] + N2[X,vib=$w]")
            )
        )
    end
end


filter!(r -> isreactive(r), rxs)

rxs = apply_tree(tr, identity.(rxs))


##

p = 1.0e5
Ngas = p / kB / 3000.0

@named sys = ReactionSystem(to_catalyst(tr, rxs...), t; combinatoric_ratelaws = false)


@named thermal = ReactionSystem(
    [
        D(Tgas) ~ sum([thermal_source_term(r)(to_catalyst.(r.formula.subs)...) for r in rxs]) / (5 / 2 * Ngas * kB),
    ];
    combinatoric_ratelaws = false
)


nsys = extend(sys, thermal)

species(nsys)

csys = complete(nsys)


##


prob = ODEProblem(
    csys,
    [
        csys.var"N2[X,vib=0]" => Ngas,
        Tgas => 3000.0,
        # csys.var"N2[X,vib=0]" => (Ngas - Nv1),
        # csys.var"N2[X,vib=1]" => Nv1,
    ], (0.0, 100.0e-3)
)

sol = solve(prob)

##

sum(sol[:, end])

Ngas


##

lines(sol.t, sol[Tgas])

##


ΔE(1) * h * 3.0e10


# lines!(
#     0 .. 10, v -> sol[2, end] * exp(-wavenumber_to_eV(energy(v)) * 1.6e-19 / kB / Tgas)
# )
