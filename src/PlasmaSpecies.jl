module PlasmaSpecies

import AbstractTrees

# include("species.jl")
# export Positive, Negative, Neutral
# export MolecularElectronicState
# export energy, charge
# export Nitrogen, DiNitrogen

include("species_node.jl")
include("species_tree.jl")

to_catalyst(_...) = error("This function needs Catalyst.jl")
parse_reaction(_...) = error("This function needs Catalyst.jl")
export to_catalyst, parse_reaction

end 
