using PlasmaSpecies
using Catalyst
using Test

@testset "Catalyst extension" begin
    @testset "to_catalyst(Species)" begin
        @test string(to_catalyst(Species("e")))  == "e(t)"
        @test string(to_catalyst(Species("N2"))) == "N2(t)"
    end

    @testset "to_catalyst(nothing) errors" begin
        @test_throws ErrorException to_catalyst(nothing)
    end

    @testset "to_catalyst(tree, (rate, formula)): no branching" begin
        tree = SpeciesTree(["N2", "e", "N2[+]"])
        result = to_catalyst(tree, (1.5e-14, ReactionFormula("N2 + e --> N2[+] + 2e")))

        @test length(result) == 1
        @test result[1] isa Reaction
        @test result[1].rate ≈ 1.5e-14
        @test length(result[1].substrates) == 2
        @test result[1].substoich == [1, 1]
        @test length(result[1].products) == 2
        @test result[1].prodstoich == [2, 1]
    end

    @testset "to_catalyst(tree, (rate, formula)): product branching" begin
        # N2[A] has one leaf, N2[B] has two leaves → two Reactions at half rate
        tree = SpeciesTree([
            "N2[A]",
            "N2[A,vib=1]",
            "N2[B]",
            "N2[B,vib=1]",
            "N2[B,vib=2]",
            "e",
        ])
        result = to_catalyst(tree, (1.0, ReactionFormula("N2[A] + e --> N2[B] + e")))

        @test length(result) == 2
        @test all(r isa Reaction for r in result)
        @test all(r.rate ≈ 0.5 for r in result)
    end

    @testset "to_catalyst(tree, rate_reactions...): variadic" begin
        tree = SpeciesTree(["N2", "e", "N2[+]"])
        rxn = ReactionFormula("N2 + e --> N2[+] + 2e")
        result = to_catalyst(tree, (1.0e-14, rxn), (2.0e-14, rxn))

        @test length(result) == 2
        rates = Set(r.rate for r in result)
        @test 1.0e-14 ∈ rates
        @test 2.0e-14 ∈ rates
    end

    @testset "to_catalyst(tree, Vector): vector dispatch" begin
        tree = SpeciesTree(["N2", "e", "N2[+]"])
        rxn = ReactionFormula("N2 + e --> N2[+] + 2e")
        result = to_catalyst(tree, [(1.0e-14, rxn), (2.0e-14, rxn)])

        @test length(result) == 2
    end

    @testset "to_catalyst(tree, PlasmaReaction): no branching" begin
        tree = SpeciesTree(["N2", "e", "N2[+]"])
        pr = PlasmaReaction(1.5e-14, ReactionFormula("N2 + e --> N2[+] + 2e"))
        result = to_catalyst(tree, pr)

        @test length(result) == 1
        @test result[1] isa Reaction
        @test result[1].rate ≈ 1.5e-14
        @test length(result[1].substrates) == 2
        @test result[1].substoich == [1, 1]
        @test length(result[1].products) == 2
        @test result[1].prodstoich == [2, 1]
    end

    @testset "to_catalyst(tree, PlasmaReaction): product branching" begin
        tree = SpeciesTree([
            "N2[A]",
            "N2[A,vib=1]",
            "N2[B]",
            "N2[B,vib=1]",
            "N2[B,vib=2]",
            "e",
        ])
        pr = PlasmaReaction(1.0, ReactionFormula("N2[A] + e --> N2[B] + e"))
        result = to_catalyst(tree, pr)

        @test length(result) == 2
        @test all(r isa Reaction for r in result)
        @test all(r.rate ≈ 0.5 for r in result)
    end

    @testset "to_catalyst(tree, PlasmaReaction...): variadic" begin
        tree = SpeciesTree(["N2", "e", "N2[+]"])
        pr1 = PlasmaReaction(1.0e-14, ReactionFormula("N2 + e --> N2[+] + 2e"))
        pr2 = PlasmaReaction(2.0e-14, ReactionFormula("N2 + e --> N2[+] + 2e"))
        result = to_catalyst(tree, pr1, pr2)

        @test length(result) == 2
        rates = Set(r.rate for r in result)
        @test 1.0e-14 ∈ rates
        @test 2.0e-14 ∈ rates
    end

    @testset "to_catalyst(tree, Vector{PlasmaReaction}): vector dispatch" begin
        tree = SpeciesTree(["N2", "e", "N2[+]"])
        prs = [
            PlasmaReaction(1.0e-14, ReactionFormula("N2 + e --> N2[+] + 2e")),
            PlasmaReaction(2.0e-14, ReactionFormula("N2 + e --> N2[+] + 2e")),
        ]
        result = to_catalyst(tree, prs)

        @test length(result) == 2
    end
end
