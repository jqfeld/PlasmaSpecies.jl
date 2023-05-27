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
    collect(Leaves(t)) .|> nodevalue .|> PlasmaSpecies.to_catalyst
end

function get_stoich(v::AbstractString)
    m = match(r"^(\d+)", v)
    return isnothing(m) ? 1 : parse(Int, m[1])
end
remove_stoich(s::AbstractString) = replace(s, r"^\d+" => s"")


function combine_same(sp, stoich)
    new_sp = empty(sp)
    new_stoich = empty(stoich)
    while !isempty(sp)
        s = pop!(sp)
        st = pop!(stoich)
        i = findfirst(==(s), new_sp)
        if isnothing(i)
            push!(new_sp, s)
            push!(new_stoich, st)
        else
            new_stoich[i] += st
        end
    end
    return new_sp, new_stoich
end


struct ReactionRecipe
    subs
    prods
    substoich
    prodstoich
    reverse::Bool
end
function Base.show(io::IO, recipe::ReactionRecipe)
    if length(recipe.subs) > 1
        reduce(function (x, y)
                if x[2] > 1
                    show(io, x[2])
                end
                show(io, x[1])
                print(io, '+')
                if y[2] > 1
                    show(io, y[2])
                end
                show(io, y[1])
            end, zip(recipe.subs, recipe.substoich))
    else
        if recipe.substoich[1] > 1
            show(io, recipe.substoich[1])
        end
        show(io, recipe.subs[1])
    end

    if recipe.reverse
        print(io, "<-->")
    else
        print(io, "-->")
    end
    
    if length(recipe.prods) > 1
        reduce(function (x, y)
                if x[2] > 1
                    show(io, x[2])
                end
                show(io, x[1])
                print(io, '+')
                if y[2] > 1
                    show(io, y[2])
                end
                show(io, y[1])
            end, zip(recipe.prods, recipe.prodstoich))
    else
        if recipe.prodstoich[1] > 1
            show(io, recipe.prodstoich[1])
        end
        show(io, recipe.prods[1])
    end
end

function PlasmaSpecies.parse_recipe(str)
    str = replace(str, r"\s" => "")
    dir = match(r"(-->|<-->|<--)", str)[1]
    reverse_reaction = dir == "<-->"
    lhs, rhs = dir == "<--" ? split(str, "<--")[end:-1:1] : split(str, "-->")
    substrates_str = split(lhs, r"(?<=[^(])\+")
    products_str = split(rhs, r"(?<=[^(])\+")
    substoich = get_stoich.(substrates_str)
    prodstoich = get_stoich.(products_str)
    subs = Species.(remove_stoich.(substrates_str))
    prods = Species.(remove_stoich.(products_str))
    subs, substoich = combine_same(subs, substoich)
    prods, prodstoich = combine_same(prods, prodstoich)
    return ReactionRecipe(
        subs,
        prods,
        substoich,
        prodstoich,
        reverse_reaction
    )
end
ReactionRecipe(str) = PlasmaSpecies.parse_recipe(str)


function PlasmaSpecies.make_reaction(t::SpeciesTree, rate_reaction::Tuple{Any,ReactionRecipe})
    rate, recipe = rate_reaction
    substrate_vectors = Iterators.product([PlasmaSpecies.to_catalyst(t[s]) for s in recipe.subs]...) .|> collect
    products_vectors = Iterators.product([PlasmaSpecies.to_catalyst(t[s]) for s in recipe.prods]...) .|> collect
    branching_factor = length(products_vectors)
    parameter_symbol = Symbol("k["*string(recipe)*"]")
[Reaction(@parameters($parameter_symbol = rate )[1]/ branching_factor, subs, prods, recipe.substoich, recipe.prodstoich) for subs in substrate_vectors for prods in products_vectors]
end

function PlasmaSpecies.make_reaction(t::SpeciesTree, rate_reactions::Tuple{Any,ReactionRecipe}...)
    out = []
    for rate_reaction in rate_reactions
        try
            append!(out, PlasmaSpecies.make_reaction(t, rate_reaction))
        catch error
            @warn "Skipping reaction because an error occured while making " rate_reaction error
        end
    end
    return out
end
PlasmaSpecies.make_reaction(t::SpeciesTree, rate_reactions::Vector{Tuple{Any,ReactionRecipe}}) = PlasmaSpecies.make_reaction(t, rate_reactions...)

end
