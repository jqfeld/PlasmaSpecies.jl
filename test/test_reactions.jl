using PlasmaSpecies
using Test

@testset "PlasmaReaction" begin
    @testset "Parsing" begin
        rxn = PlasmaReaction("N2 + e --> N2(+) + 2e")
        @test !rxn.reverse

        substrate_counts = Dict(string.(rxn.subs) .=> rxn.substoich)
        product_counts = Dict(string.(rxn.prods) .=> rxn.prodstoich)

        @test substrate_counts == Dict("N2" => 1, "e" => 1)
        @test product_counts == Dict("N2(+)" => 1, "e" => 2)

        reversible = PlasmaReaction("N2 + e <--> N2(+) + 2e")
        @test reversible.reverse
        @test Set(string.(reversible.subs)) == Set(["N2", "e"])
        @test Set(string.(reversible.prods)) == Set(["N2(+)", "e"])

        backward = PlasmaReaction("N2 <-- N2(+) + e")
        @test !backward.reverse
        @test Dict(string.(backward.subs) .=> backward.substoich) == Dict("N2(+)" => 1, "e" => 1)
        @test Dict(string.(backward.prods) .=> backward.prodstoich) == Dict("N2" => 1)
    end

    @testset "Display" begin
        rxn = PlasmaReaction("N2 + e --> N2(+) + e")
        @test string(rxn) == "e+N2-->e+N2(+)"
        @test sprint(show, rxn) == "e+N2-->e+N2(+)"
    end

    @testset "Macro" begin
        @test p"N2" == Species("N2")
        macro_rxn = p"N2 + e --> N2(+) + 2e"
        @test macro_rxn isa PlasmaReaction
        @test macro_rxn.prodstoich == PlasmaReaction("N2 + e --> N2(+) + 2e").prodstoich
    end

    @testset "Mass balance" begin
        balanced = PlasmaReaction("N2 + e --> N2(+) + 2e")
        @test PlasmaSpecies.ismassbalanced(balanced)

        unbalanced = PlasmaReaction("N2 + e --> N2(+) + e")
        @test !PlasmaSpecies.ismassbalanced(unbalanced)
    end

    @testset "Species tree expansion" begin
        tree = SpeciesTree([
            "N2(A)",
            "N2(A,v=1)",
            "N2(B)",
            "N2(B,v=1)",
            "N2(B,v=2)",
            "e",
        ])

        rxn = PlasmaReaction("N2(A) + e --> N2(B) + e")
        expanded = apply_tree(tree, rxn)

        @test length(expanded) == 2
        @test all(first.(expanded) .== 0.5)

        expanded_rxns = last.(expanded)
        @test all(string.(rxn.subs) == ["e", "N2(A,v=1)"] for rxn in expanded_rxns)

        product_signatures = Set(Tuple(string.(rxn.prods)) for rxn in expanded_rxns)
        @test product_signatures == Set([
            ("e", "N2(B,v=1)"),
            ("e", "N2(B,v=2)"),
        ])
    end
end
