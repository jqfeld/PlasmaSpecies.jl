module PlasmaSpecies

import AbstractTrees
include("species.jl")
include("species_tree.jl")


to_catalyst(_) = error("This function need Catalyst.jl")
export to_catalyst
end
