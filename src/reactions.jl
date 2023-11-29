
struct PlasmaReaction
    subs
    prods
    substoich
    prodstoich
    reverse::Bool
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


function Base.show(io::IO, recipe::PlasmaReaction)
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

function sort_by_species(sp, stoich)
    p = sortperm(sp, rev=true)
    return sp[p], stoich[p]
end

function parse_reaction(str)
    str = replace(str, r"\s" => "")
    dir = match(r"(-->|<-->|<--)", str)[1]
    reverse_reaction = dir == "<-->"
    lhs, rhs = dir == "<--" ? split(str, "<--")[end:-1:1] : split(str, "-->")
    substrates_str = split(lhs, r"(?<!\(|,)\+(?!,|\)|\s\)|\s,)")
    products_str = split(rhs, r"(?<!\(|,)\+(?!,|\)|\s\)|\s,)")
    substoich = get_stoich.(substrates_str)
    prodstoich = get_stoich.(products_str)
    subs = Species.(remove_stoich.(substrates_str))
    prods = Species.(remove_stoich.(products_str))
    subs, substoich = sort_by_species(combine_same(subs, substoich)...)
    prods, prodstoich = sort_by_species(combine_same(prods, prodstoich)...)
    return PlasmaReaction(
        subs,
        prods,
        substoich,
        prodstoich,
        reverse_reaction
    )
end
PlasmaReaction(str) = parse_reaction(str)

macro p_str(s)
    if contains(s, r"(-->|<-->|<--)")
        return parse_reaction(s)
    else
        return Species(s)
    end

end

"""
    ```julia 
    apply_tree(t::SpeciesTree, reactions::PlasmaReaction)
    ```

Apply the species tree to a reaction and return an array with tuples containing 
the scaling factor and the new reaction. 
This means for each participating species it is checked, if it is a leaf of the tree. 
If not, a new reaction is created for each descendent species. 
If the non-leaf species is a product, the branching factor is one over the number of 
descendants (effectively assuming that the total reaction rate is distributed equally
over all possible products).

"""
function apply_tree(t::SpeciesTree, reaction::PlasmaReaction)
    substrate_vectors = Iterators.product([leaves(t[s]) for s in reaction.subs]...) .|> collect
    products_vectors = Iterators.product([leaves(t[s]) for s in reaction.prods]...) .|> collect
    branching_factor = 1 / length(products_vectors)
    [(branching_factor, PlasmaReaction(subs, prods, reaction.substoich, reaction.prodstoich, reaction.reverse)) for subs in substrate_vectors for prods in products_vectors]
end

function apply_tree(t::SpeciesTree, reactions::PlasmaReaction...)
    out = []
    for rate_reaction in reactions
        try
            append!(out, apply_tree(t, rate_reaction))
        catch error
            @warn "Skipping reaction because an error occured while making " rate_reaction error
        end
    end
    return out
end
apply_tree(t::SpeciesTree, reactions::Vector{PlasmaReaction}) = PlasmaSpecies.apply_tree(t, reactions...)


ismassbalanced(r::PlasmaReaction) = sum(mass.(r.prods) .* r.prodstoich) ≈ sum(mass.(r.subs) .* r.substoich)
ischargebalanced(r::PlasmaReaction) = sum(to_value.(charge.(r.prods)) .* r.prodstoich) == sum(to_value.(charge.(r.subs)) .* r.substoich)
