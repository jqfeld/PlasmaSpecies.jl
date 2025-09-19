using SafeTestsets

@time begin
    # Write your tests here.
    @time @safetestset "Species" begin include("test_species.jl") end
    @time @safetestset "SpeciesNode" begin include("test_species_node.jl") end
    @time @safetestset "SpeciesTree" begin include("test_species_tree.jl") end
    @time @safetestset "mass" begin include("test_mass.jl") end
    @time @safetestset "reactions" begin include("test_reactions.jl") end
end
