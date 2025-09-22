using PlasmaSpecies
using Test


@test mass(Species("N2")) == mass(DiNitrogen())
@test mass(Species("N2(X,v=0)")) == mass(Species("N2"))
@test mass(Species("N2(X,v=0)")) == mass(p"N2")
@test mass(Species("N2(-)")) ≈ mass(Species("N2")) + mass(Species("e"))
@test mass(Species("N2(+)")) ≈ mass(Species("N2")) - mass(Species("e"))

@test mass(Species("O2")) == mass(DiOxygen())
@test mass(Species("O2(X,v=0)")) == mass(Species("O2"))
@test mass(Species("O2(X,v=0)")) == mass(p"O2")
@test mass(Species("O2(-)")) ≈ mass(Species("O2")) + mass(Species("e"))
@test mass(Species("O2(+)")) ≈ mass(Species("O2")) - mass(Species("e"))

