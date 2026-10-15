using PlasmaSpecies
using Test
using AbstractTrees

# @test unique( x-> nodevalue(x), children(combine!(SpeciesTree("N2[X,vib=0,rot=0]"), SpeciesTree("N2[X,vib=0,rot=1]"))))

@testset "SpeciesTree" begin
    species = [
        Species("N2[A]"),
        Species("N2[A,vib=1]"),
        Species("N2[B]"),
        Species("N2[B,vib=1]"),
        Species("e"),
    ]

    @testset "Vector{Species} matches Vector{String}" begin
        tree_from_species = SpeciesTree(species)
        tree_from_strings = SpeciesTree(string.(species))
        @test leaves(tree_from_species) == leaves(tree_from_strings)
    end

    @testset "variadic Species... matches variadic String..." begin
        tree_from_species = SpeciesTree(species...)
        tree_from_strings = SpeciesTree(string.(species)...)
        @test leaves(tree_from_species) == leaves(tree_from_strings)
    end

    @testset "leaves are correct" begin
        tree = SpeciesTree(species)
        ls = leaves(tree)
        @test Species("N2[A,vib=1]") ∈ ls
        @test Species("N2[B,vib=1]") ∈ ls
        @test Species("e") ∈ ls
        @test Species("N2[A]") ∉ ls
        @test Species("N2[B]") ∉ ls
    end
end

@testset "empty and incremental trees" begin
    for t in (SpeciesTree(), SpeciesTree(String[]), SpeciesTree(Species[]))
        @test isnothing(t.x)
        @test isempty(leaves(t))
    end
    t = SpeciesTree()
    push!(t, "N2[X,vib=0]")
    push!(t, Species("N2[X,vib=1]"))
    @test string.(leaves(t)) == ["N2[X,vib=0]", "N2[X,vib=1]"]
    # SubString species labels
    @test leaves(SpeciesTree(split("e N2[+]"))) == Species.(["e", "N2[+]"])
    # an empty tree survives the energy pass
    @test isempty(leaves(apply_energy!(SpeciesTree(), _ -> 0.0)))
end
