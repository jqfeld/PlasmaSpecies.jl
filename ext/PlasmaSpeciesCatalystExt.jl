module PlasmaSpeciesCatalystExt

using AbstractTrees
using PlasmaSpecies
using Catalyst

function PlasmaSpecies.to_catalyst(sp::Species)
    # Catalyst's shared independent variable, so these species match the `t`
    # of any ReactionSystem they are put into.
    t = default_t()
    symbol = Symbol(string(sp))
    return (@species $symbol(t) = 0.0)[1]
end

function PlasmaSpecies.to_catalyst(::Nothing)
    error("Species not found.")
end

function PlasmaSpecies.to_catalyst(t::SpeciesTree)
    leaves(t) .|> PlasmaSpecies.to_catalyst |> sum
end

function PlasmaSpecies.to_catalyst(t::SpeciesTree, rate_reaction::Tuple{Any,ReactionFormula})
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

# Batches: one typed vector, and an error in any reaction aborts the batch
# (logged with the reaction, then rethrown) rather than being skipped.
function PlasmaSpecies.to_catalyst(t::SpeciesTree, rate_reactions::Tuple{Any,ReactionFormula}...)
    out = Reaction[]
    for rate_reaction in rate_reactions
        append!(out, PlasmaSpecies._expanding(() -> PlasmaSpecies.to_catalyst(t, rate_reaction), rate_reaction))
    end
    return out
end
PlasmaSpecies.to_catalyst(t::SpeciesTree, rate_reactions::AbstractVector{<:Tuple{Any,ReactionFormula}}) =
    PlasmaSpecies.to_catalyst(t, rate_reactions...)

function PlasmaSpecies.to_catalyst(t::SpeciesTree, pr::PlasmaReaction)
    [
        Reaction(
            r.rate,
            PlasmaSpecies.to_catalyst.(r.formula.subs),
            PlasmaSpecies.to_catalyst.(r.formula.prods),
            r.formula.substoich,
            r.formula.prodstoich
        ) for r in apply_tree(t, pr)
    ]
end

function PlasmaSpecies.to_catalyst(t::SpeciesTree, prs::PlasmaReaction...)
    out = Reaction[]
    for pr in prs
        append!(out, PlasmaSpecies._expanding(() -> PlasmaSpecies.to_catalyst(t, pr), pr))
    end
    return out
end
PlasmaSpecies.to_catalyst(t::SpeciesTree, prs::AbstractVector{<:PlasmaReaction}) =
    PlasmaSpecies.to_catalyst(t, prs...)

end
