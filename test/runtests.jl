using SafeTestsets

@time begin
    # Write your tests here.
    @time @safetestset "Species" begin include("test_species.jl") end 
    @time @safetestset "SpeciesTree" begin include("test_species_tree.jl") end 
end
