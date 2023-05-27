module PlasmaSpecies

import AbstractTrees

# include("species.jl")
export Positive, Negative, Neutral, Charge
# export MolecularElectronicState
# export energy, charge
# export Nitrogen, DiNitrogen

include("species_node.jl")
include("species_tree.jl")

to_catalyst(_...) = error("This function needs Catalyst.jl")
make_reaction(_...) = error("This function needs Catalyst.jl")
parse_recipe(_...) = error("This function needs Catalyst.jl")

macro p_str(s)
    parse_recipe(s)
end
export @p_str, parse_recipe

export to_catalyst, make_reaction

end 
