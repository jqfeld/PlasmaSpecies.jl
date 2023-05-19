using PlasmaSpecies
using Test


@test unique( x-> nodevalue(x), children(combine(SpeciesTree("N2(X,v=0,J=0)"), SpeciesTree("N2(X,v=0,J=1)"))))
