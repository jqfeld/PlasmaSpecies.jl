using PlasmaSpecies
using Test

@test get_parent_species(SpeciesNode("N2(X,v=0)")) == SpeciesNode("N2(X)")
@test get_parent_species(SpeciesNode("N2(+,X)")) == SpeciesNode("N2(+)")
@test get_parent_species(SpeciesNode("N2(+)")) !== SpeciesNode("N2")
@test get_parent_species(SpeciesNode("N2(X,v=2,J=10)")) == SpeciesNode("N2(X,v=2)")
@test get_parent_species(SpeciesNode("N2(X,v=2,J=10)")) !== SpeciesNode("N2(X,v=1)")
@test is_parent_species(get_parent_species(SpeciesNode("N2(X,v=2,J=10)")), SpeciesNode("N2(X)"))
@test !is_parent_species(get_parent_species(SpeciesNode("N2(X,v=2,J=10)")), SpeciesNode("N2(A)"))
@test !is_parent_species(get_parent_species(SpeciesNode("N2(X,v=2,J=10)")), SpeciesNode("N2(a'')"))
@test !is_parent_species(get_parent_species(SpeciesNode("N2(X,v=2,J=10)")), SpeciesNode("N(4S)"))
