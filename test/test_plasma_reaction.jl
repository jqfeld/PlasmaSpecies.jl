using PlasmaSpecies
using Test

@testset "PlasmaReaction" begin
    @testset "Construction" begin
        formula = ReactionFormula("N2 + e --> N2[+] + 2e")
        pr = PlasmaReaction(1.5e-14, formula)
        @test pr.rate == 1.5e-14
        @test pr.formula === formula
        # the formula compares by value, so an equivalent one built elsewhere matches
        @test pr.formula == ReactionFormula("e + N2 --> 2e + N2[+]")
    end

    @testset "apply_tree: no branching" begin
        # All species are leaves with no children, so branching_factor = 1
        tree = SpeciesTree(["N2", "e", "N2[+]"])
        pr = PlasmaReaction(1.0e-14, ReactionFormula("N2 + e --> N2[+] + 2e"))
        expanded = apply_tree(tree, pr)

        @test length(expanded) == 1
        @test expanded[1] isa PlasmaReaction
        @test expanded[1].rate ≈ 1.0e-14
        @test string.(expanded[1].formula.subs) == ["e", "N2"]
        @test string.(expanded[1].formula.prods) == ["e", "N2[+]"]
    end

    @testset "apply_tree: product branching scales rate" begin
        # N2[A] has one leaf, N2[B] has two leaves → branching_factor = 0.5
        tree = SpeciesTree([
            "N2[A]",
            "N2[A,vib=1]",
            "N2[B]",
            "N2[B,vib=1]",
            "N2[B,vib=2]",
            "e",
        ])
        pr = PlasmaReaction(1.0, ReactionFormula("N2[A] + e --> N2[B] + e"))
        expanded = apply_tree(tree, pr)

        @test length(expanded) == 2
        @test all(r isa PlasmaReaction for r in expanded)
        @test all(r.rate ≈ 0.5 for r in expanded)

        product_signatures = Set(Tuple(string.(r.formula.prods)) for r in expanded)
        @test product_signatures == Set([
            ("e", "N2[B,vib=1]"),
            ("e", "N2[B,vib=2]"),
        ])
    end

    @testset "energy: known energies" begin
        sp1 = Species(gas=Gas("N2"), energy=2.0)
        sp2 = Species(gas=Gas("N2"), electronic_state=StringElectronicState("A"), energy=5.0)
        pr = PlasmaReaction(1.0, ReactionFormula([sp1], [sp2], [1], [1], false))
        @test reaction_energy(pr) ≈ 3.0
    end

    @testset "energy: stoichiometry" begin
        sp1 = Species(gas=Gas("N2"), energy=1.0)
        sp2 = Species(gas=Gas("N2"), electronic_state=StringElectronicState("A"), energy=4.0)
        # 2*sp1 --> sp2: energy = 4.0 - 2*1.0 = 2.0
        pr = PlasmaReaction(1.0, ReactionFormula([sp1], [sp2], [2], [1], false))
        @test reaction_energy(pr) ≈ 2.0
    end

    @testset "energy: unknown energy returns nothing" begin
        sp1 = Species(gas=Gas("N2"))
        sp2 = Species(gas=Gas("N2"), electronic_state=StringElectronicState("A"), energy=5.0)
        pr_unknown_sub = PlasmaReaction(1.0, ReactionFormula([sp1], [sp2], [1], [1], false))
        @test isnothing(reaction_energy(pr_unknown_sub))

        sp3 = Species(gas=Gas("N2"), energy=2.0)
        pr_unknown_prod = PlasmaReaction(1.0, ReactionFormula([sp3], [sp1], [1], [1], false))
        @test isnothing(reaction_energy(pr_unknown_prod))
    end

    @testset "apply_tree: Vector{PlasmaReaction}" begin
        tree = SpeciesTree(["N2", "e", "N2[+]"])
        prs = [
            PlasmaReaction(1.0e-14, ReactionFormula("N2 + e --> N2[+] + 2e")),
            PlasmaReaction(2.0e-14, ReactionFormula("N2 + e --> N2 + e")),
        ]
        expanded = apply_tree(tree, prs)
        @test length(expanded) == 2
        @test expanded[1].rate ≈ 1.0e-14
        @test expanded[2].rate ≈ 2.0e-14
    end

    @testset "isreactive: different LHS and RHS" begin
        @test isreactive(PlasmaReaction(1.0, ReactionFormula("N2 + e --> N2[+] + 2e")))
        @test isreactive(ReactionFormula("N2 + e --> N2[+] + 2e"))
    end

    @testset "isreactive: identical LHS and RHS" begin
        @test !isreactive(PlasmaReaction(1.0, ReactionFormula("N2 + e --> e + N2")))
        @test !isreactive(ReactionFormula("N2 + e --> e + N2"))
    end

    @testset "isreactive: same species different stoichiometry" begin
        sp = Species(gas=Gas("N2"), energy=1.0)
        @test isreactive(ReactionFormula([sp], [sp], [1], [2], false))
    end

    lv(s, E) = with_fields(Species(s), energy=E)
    el = Species("e")
    n2 = lv("N2[X]", 0.0); n2v = lv("N2[X,vib=1]", 0.29)
    he = lv("He[X]", 0.0); hes = lv("He[*]", 19.82); hep = lv("He[+]", 24.58)
    he2p = lv("He2[+]", 23.03); he2s = lv("He2[*]", 18.4); n2p = lv("N2[+]", 15.58)
    rx(subs, prods, ss, ps; kw...) = PlasmaReaction(1.0, ReactionFormula(subs, prods, ss, ps, false); kw...)

    @testset "energy_partition: the rule" begin
        vt = rx([n2v, n2], [n2], [1, 1], [2])                       # VT: two heavy products, the gas
        @test energy_partition(vt).gas ≈ 0.29 && energy_partition(vt).electrons == 0 && energy_partition(vt).radiation == 0
        pen = rx([hes, n2], [he, n2p, el], [1, 1], [1, 1, 1])      # Penning: an electron product, the electrons
        @test energy_partition(pen).electrons ≈ 19.82 - 15.58 && energy_partition(pen).gas == 0
        ion = rx([el, hes], [el, hep], [1, 1], [2, 1])             # stepwise ionisation: the electrons pay
        @test energy_partition(ion).electrons ≈ 19.82 - 24.58
        rr = rx([el, hep], [hes], [1, 1], [1])                     # two-body recombination: one product, radiation
        p = energy_partition(rr; consumed=2.0)
        @test p.radiation ≈ 24.58 - 19.82 + 2.0 && p.electrons == 0 && p.gas == 0   # the electron's kinetic energy too
        dr = rx([el, he2p], [hes, he], [1, 1], [1, 1])             # dissociative recombination: the fragments, the gas
        @test energy_partition(dr; consumed=1.0).gas ≈ 23.03 - 19.82 + 1.0
        tb = rx([el, hep], [el, hes], [2, 1], [1, 1])              # three-body on an electron: the electrons
        @test energy_partition(tb; consumed=1.0).electrons ≈ 24.58 - 19.82 + 1.0
        @test energy_partition(pen; consumed=5.0) == energy_partition(pen)   # it consumes no electron
        @test energy_partition(rx([lv("N2[A]", nothing), n2], [n2], [1, 1], [2])) === nothing   # a level unknown
        @test reaction_energy(pen) ≈ 15.58 - 19.82                 # an electron's level: 0
    end

    @testset "energy_partition: the reaction's own shares" begin
        decay = rx([he2s], [he], [1], [2]; partition=(radiation=1,))   # an excimer's decay radiates
        @test energy_partition(decay) == (electrons=0.0, gas=0.0, radiation=18.4)
        @test energy_partition(rx([he2s], [he], [1], [2])).gas ≈ 18.4      # the rule: two bodies, the gas
        half = rx([n2v, n2], [n2], [1, 1], [2]; partition=(gas=0.5, radiation=0.5))
        @test energy_partition(half).gas ≈ 0.145 && energy_partition(half).radiation ≈ 0.145
        @test_throws ArgumentError rx([he2s], [he], [1], [2]; partition=(gas=0.5,))      # not summing to 1
        @test_throws ArgumentError rx([he2s], [he], [1], [2]; partition=(heat=1.0,))     # no such channel
        @test PlasmaReaction(1.0, decay.formula).partition === nothing
        # apply_tree splits the reaction into states: each child keeps the shares
        tree = SpeciesTree([he2s, he])
        @test all(c -> c.partition == decay.partition, apply_tree(tree, decay))
    end

    @testset "energy_source_terms and thermal_source_term: the channels' power" begin
        pen = PlasmaReaction(2.0, ReactionFormula([hes, n2], [he, n2p, el], [1, 1], [1, 1, 1], false))
        dens = Dict(string(hes) => 3.0, string(n2) => 4.0)
        ns = [dens[string(sp)] for sp in pen.formula.subs]
        q = energy_source_terms(pen)(ns...)
        @test q.electrons ≈ 2.0 * 12.0 * (19.82 - 15.58) && q.gas == 0 && q.radiation == 0
        @test thermal_source_term(pen)(ns...) == 0                # the electron takes it, not the gas (it used to)
        vt = PlasmaReaction(Tg -> 1e-3 * Tg, ReactionFormula([n2v, n2], [n2], [1, 1], [2], false))
        f = thermal_source_term(vt)
        @test f(2.0, 5.0; rate_args=(300.0,)) ≈ 0.3 * 10.0 * 0.29    # a callable rate, its arguments
        @test_throws MethodError f(2.0, 5.0)                          # it needs them
        @test_throws ArgumentError f(2.0; rate_args=(300.0,))         # the substrates' count
        dr = PlasmaReaction(1.0, ReactionFormula([el, he2p], [hes, he], [1, 1], [1, 1], false))
        @test thermal_source_term(dr)(1.0, 1.0; consumed=1.0) ≈ 23.03 - 19.82 + 1.0
        @test_throws ErrorException thermal_source_term(PlasmaReaction(1.0, ReactionFormula([lv("N2[A]", nothing), n2], [n2], [1, 1], [2], false)))
    end

    @testset "thermal_source_term: unknown energy throws" begin
        sp1 = Species(gas=Gas("N2"))
        sp2 = Species(gas=Gas("N2"), electronic_state=StringElectronicState("A"), energy=5.0)
        pr = PlasmaReaction(1.0, ReactionFormula([sp1], [sp2], [1], [1], false))
        @test_throws ErrorException thermal_source_term(pr)
    end

    @testset "apply_tree: rate is distributed across branches" begin
        # N2[A,vib=1] is a leaf (substrate, no branching); N2[B] has two leaves (product, branching_factor = 0.5)
        # → 2 expanded reactions each with rate 6.0 * 0.5 = 3.0, sum = 6.0
        tree = SpeciesTree([
            "N2[A,vib=1]",
            "N2[B]",
            "N2[B,vib=1]",
            "N2[B,vib=2]",
            "e",
        ])
        pr = PlasmaReaction(6.0, ReactionFormula("N2[A,vib=1] + e --> N2[B] + e"))
        expanded = apply_tree(tree, pr)

        @test length(expanded) == 2
        @test all(r.rate ≈ 3.0 for r in expanded)
        @test sum(r.rate for r in expanded) ≈ 6.0
    end
end
