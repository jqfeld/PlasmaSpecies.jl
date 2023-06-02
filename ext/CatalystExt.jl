module CatalystExt

using AbstractTrees
using PlasmaSpecies
using Catalyst

function PlasmaSpecies.to_catalyst(sp::Species)
    @variables t
    symbol = Symbol(string(sp))
    return (@species $symbol(t) = 0.0)[1]
end

function PlasmaSpecies.to_catalyst(::Nothing)
    error("Species not found.")
end

function PlasmaSpecies.to_catalyst(t::SpeciesTree)
    leaves(t) .|> PlasmaSpecies.to_catalyst |> sum
end

function PlasmaSpecies.to_catalyst(t::SpeciesTree, rate_reaction::Tuple{Any,PlasmaReaction})
    rate, reaction = rate_reaction
    scaled_reactions = apply_tree(t, reaction)
    # parameter_symbol = Symbol("k["*string(reaction)*"]")
    return [
        Reaction(
            rate * branching_factor, 
            PlasmaSpecies.to_catalyst.(r.subs), 
            PlasmaSpecies.to_catalyst.(r.prods), 
            r.substoich, r.prodstoich
        ) for (branching_factor, r) in scaled_reactions
    ]
end

function PlasmaSpecies.to_catalyst(t::SpeciesTree, rate_reactions::Tuple{Any,PlasmaReaction}...)
    out = []
    for rate_reaction in rate_reactions
        try
            append!(out, PlasmaSpecies.to_catalyst(t, rate_reaction))
        catch error
            @warn "Skipping reaction because an error occured while making " rate_reaction error
        end
    end
    return out
end
PlasmaSpecies.to_catalyst(t::SpeciesTree, rate_reactions::Vector{Tuple{Any,PlasmaReaction}}) = PlasmaSpecies.to_catalyst(t, rate_reactions...)

end
