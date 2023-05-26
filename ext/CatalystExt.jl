module CatalystExt

using AbstractTrees
using PlasmaSpecies
using Catalyst

function PlasmaSpecies.to_catalyst(sp::SpeciesNode)
    @variables t
    symbol = Symbol(string(sp))
    return (@species $symbol(t) = 0.0)[1]
end

function PlasmaSpecies.to_catalyst(::Nothing)
    error("Species not found.")
end

function PlasmaSpecies.to_catalyst(t::SpeciesTree)
    collect(Leaves(t)) .|> nodevalue .|> PlasmaSpecies.to_catalyst
end

function get_stoich(v::AbstractString)
    m = match(r"^(\d+)", v)
    return isnothing(m) ? 1 : parse(Int, m[1])
end
remove_stoich(s::AbstractString) = replace(s, r"^\d+" => s"")

function PlasmaSpecies.parse_reaction(t::SpeciesTree, rate_reaction::Tuple{Any,String})
    rate, str = rate_reaction
    str = replace(str, r"\s" => "")
    dir = match(r"(-->|<-->|<--)", str)[1]
    reverse_reaction = dir == "<-->"
    lhs, rhs = dir == "<--" ? split(str, "<--")[end:-1:1] : split(str, "-->")
    substrates_str = split(lhs, r"(?<=[^(])\+")
    products_str = split(rhs, r"(?<=[^(])\+")
    substoich = get_stoich.(substrates_str)
    prodstoich = get_stoich.(products_str)
    substrates_str = remove_stoich.(substrates_str)
    products_str = remove_stoich.(products_str)

    substrate_vectors = Iterators.product([PlasmaSpecies.to_catalyst(t[s]) for s in substrates_str]...) .|> collect
    products_vectors = Iterators.product([PlasmaSpecies.to_catalyst(t[s]) for s in products_str]...) .|> collect
    branching_factor = length(products_vectors)

    # [println(rate /branching_factor, subs, prods, substoich, prodstoich) for subs in substrate_vectors, prods in products_vectors]
    [Reaction(rate / branching_factor, subs, prods, substoich, prodstoich) for subs in substrate_vectors, prods in products_vectors]
end

function PlasmaSpecies.parse_reaction(t::SpeciesTree, rate_reactions::Tuple{Any,String}...)
    out = []
    for rate_reaction in rate_reactions
        append!(out, PlasmaSpecies.parse_reaction(t, rate_reaction))
    end
    return out
end

end
