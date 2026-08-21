using PlasmaSpecies
using Test

@testset "ReactionFormula" begin
    @testset "Parsing" begin
        rxn = ReactionFormula("N2 + e --> N2[+] + 2e")
        @test !rxn.reverse

        substrate_counts = Dict(string.(rxn.subs) .=> rxn.substoich)
        product_counts = Dict(string.(rxn.prods) .=> rxn.prodstoich)

        @test substrate_counts == Dict("N2" => 1, "e" => 1)
        @test product_counts == Dict("N2[+]" => 1, "e" => 2)

        reversible = ReactionFormula("N2 + e <--> N2[+] + 2e")
        @test reversible.reverse
        @test Set(string.(reversible.subs)) == Set(["N2", "e"])
        @test Set(string.(reversible.prods)) == Set(["N2[+]", "e"])

        backward = ReactionFormula("N2 <-- N2[+] + e")
        @test !backward.reverse
        @test Dict(string.(backward.subs) .=> backward.substoich) == Dict("N2[+]" => 1, "e" => 1)
        @test Dict(string.(backward.prods) .=> backward.prodstoich) == Dict("N2" => 1)
    end

    @testset "Display" begin
        rxn = ReactionFormula("N2 + e --> N2[+] + e")
        @test string(rxn) == "e+N2-->e+N2[+]"
        @test sprint(show, rxn) == "e+N2-->e+N2[+]"
    end

    @testset "Macro" begin
        @test p"N2" == Species("N2")
        macro_rxn = p"N2 + e --> N2[+] + 2e"
        @test macro_rxn isa ReactionFormula
        @test macro_rxn.prodstoich == ReactionFormula("N2 + e --> N2[+] + 2e").prodstoich
    end

    @testset "Mass balance" begin
        balanced = ReactionFormula("N2 + e --> N2[+] + 2e")
        @test PlasmaSpecies.ismassbalanced(balanced)

        unbalanced = ReactionFormula("N2 + e --> N2[+] + e")
        @test !PlasmaSpecies.ismassbalanced(unbalanced)
    end

    @testset "Species tree expansion" begin
        tree = SpeciesTree([
            "N2[A]",
            "N2[A,vib=1]",
            "N2[B]",
            "N2[B,vib=1]",
            "N2[B,vib=2]",
            "e",
        ])

        rxn = ReactionFormula("N2[A] + e --> N2[B] + e")
        expanded = apply_tree(tree, rxn)

        @test length(expanded) == 2
        @test all(first.(expanded) .== 0.5)

        expanded_rxns = last.(expanded)
        @test all(string.(r.subs) == ["e", "N2[A,vib=1]"] for r in expanded_rxns)

        product_signatures = Set(Tuple(string.(rxn.prods)) for rxn in expanded_rxns)
        @test product_signatures == Set([
            ("e", "N2[B,vib=1]"),
            ("e", "N2[B,vib=2]"),
        ])
    end

    @testset "Species tree expansion: vibrational range" begin
        # vib=0-4 in the reaction should match every vib level 0..4 actually
        # present under N2[X], not just an opaque "0-4" label — vib=5 must be excluded.
        tree = SpeciesTree([
            "N2[X,vib=0]", "N2[X,vib=1]", "N2[X,vib=2]", "N2[X,vib=3]", "N2[X,vib=4]",
            "N2[X,vib=5]",
            "N2[A]",
            "e",
        ])

        rxn = ReactionFormula("N2[X,vib=0-4] + e --> N2[A] + e")
        expanded = apply_tree(tree, rxn)

        @test length(expanded) == 5
        substrate_signatures = Set(Tuple(string.(r.subs)) for (_, r) in expanded)
        @test substrate_signatures == Set([
            ("e", "N2[X,vib=$v]") for v in 0:4
        ])
        @test all(string.(r.prods) == ["e", "N2[A]"] for (_, r) in expanded)
    end
end

@testset "arrow forms" begin
    fwd = ReactionFormula("e + N2 --> 2e + N2[+]")
    rev = ReactionFormula("e + N2 <--> 2e + N2[+]")
    back = ReactionFormula("2e + N2[+] <-- e + N2")

    @test string.(fwd.subs) == ["e", "N2"]
    @test string.(fwd.prods) == ["e", "N2[+]"]
    @test !fwd.reverse

    # The reversible arrow must not leave a stray '<' on the left-hand side.
    @test string.(rev.subs) == ["e", "N2"]
    @test rev.reverse

    @test string.(back.subs) == ["e", "N2"]
    @test string.(back.prods) == ["e", "N2[+]"]
    @test !back.reverse
end

@testset "apply_tree robustness" begin
    tree = SpeciesTree(["N2[X,vib=0]", "N2[X,vib=1]", "N2[X,vib=2]", "e"])

    # A species missing from the tree drops the reaction, with a warning rather
    # than an empty result nobody notices.
    @test_logs (:warn,) match_mode = :any begin
        @test isempty(apply_tree(tree, ReactionFormula("e + Ar --> e + Ar[+]")))
    end

    # Branching factors over the products of one substrate combination sum to 1.
    expanded = apply_tree(tree, ReactionFormula("N2[X,vib=0] --> N2[X]"))
    @test length(expanded) == 3
    @test sum(first.(expanded)) ≈ 1.0

    # A callable rate is scaled, not multiplied.
    callable = PlasmaReaction(Te -> 1.0e-15 * Te, ReactionFormula("N2[X,vib=0] --> N2[X]"))
    scaled = apply_tree(tree, callable)
    @test length(scaled) == 3
    @test scaled[1].rate(2.0) ≈ 1.0e-15 * 2.0 / 3

    numeric = PlasmaReaction(1.0e-15, ReactionFormula("N2[X,vib=0] --> N2[X]"))
    @test apply_tree(tree, numeric)[1].rate ≈ 1.0e-15 / 3
end

@testset "tree indexing errors" begin
    tree = SpeciesTree(["N2[X,vib=0]", "e"])
    @test tree[Species("Ar")] === nothing
    @test_throws KeyError tree[Species("Ar")] = SpeciesTree("Ar")
end
