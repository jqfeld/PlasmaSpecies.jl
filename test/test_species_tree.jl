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
