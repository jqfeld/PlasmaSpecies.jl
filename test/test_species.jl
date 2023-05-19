using PlasmaSpecies
using Test

@test get_parent_species(Species("N2(X,v=0)")) == Species("N2(X)")
@test get_parent_species(Species("N2(+)")) == Species("N2(+,X)")
@test get_parent_species(Species("N2")) !== Species("N2(+)")
@test get_parent_species(Species("N2(X,v=2,J=10)")) == Species("N2(X,v=2)")
@test get_parent_species(Species("N2(X,v=2,J=10)")) !== Species("N2(X,v=1)")
@test is_parent_species(get_parent_species(Species("N2(X,v=2,J=10)")), Species("N2(X)"))
@test !is_parent_species(get_parent_species(Species("N2(X,v=2,J=10)")), Species("N2(A)"))
@test !is_parent_species(get_parent_species(Species("N2(X,v=2,J=10)")), Species("N2(a'')"))
@test !is_parent_species(get_parent_species(Species("N2(X,v=2,J=10)")), Species("N(4S)"))
