using PlasmaSpecies
using Test

@testset "PlasmaReaction" begin
    @testset "Construction" begin
        formula = ReactionFormula("N2 + e --> N2[+] + 2e")
        pr = PlasmaReaction(1.5e-14, formula)
        @test pr.rate == 1.5e-14
        @test pr.formula === formula
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
        sp1 = Species(gas=DiNitrogen(), energy=2.0)
        sp2 = Species(gas=DiNitrogen(), electronic_state=StringElectronicState("A"), energy=5.0)
        pr = PlasmaReaction(1.0, ReactionFormula([sp1], [sp2], [1], [1], false))
        @test reaction_energy(pr) ≈ 3.0
    end

    @testset "energy: stoichiometry" begin
        sp1 = Species(gas=DiNitrogen(), energy=1.0)
        sp2 = Species(gas=DiNitrogen(), electronic_state=StringElectronicState("A"), energy=4.0)
        # 2*sp1 --> sp2: energy = 4.0 - 2*1.0 = 2.0
        pr = PlasmaReaction(1.0, ReactionFormula([sp1], [sp2], [2], [1], false))
        @test reaction_energy(pr) ≈ 2.0
    end

    @testset "energy: unknown energy returns nothing" begin
        sp1 = Species(gas=DiNitrogen())
        sp2 = Species(gas=DiNitrogen(), electronic_state=StringElectronicState("A"), energy=5.0)
        pr_unknown_sub = PlasmaReaction(1.0, ReactionFormula([sp1], [sp2], [1], [1], false))
        @test isnothing(reaction_energy(pr_unknown_sub))

        sp3 = Species(gas=DiNitrogen(), energy=2.0)
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
