using SafeTestsets

@time begin
    # Write your tests here.
    @time @safetestset "SpeciesNode" begin include("test_species_node.jl") end 
    @time @safetestset "SpeciesTree" begin include("test_species_tree.jl") end 
end
