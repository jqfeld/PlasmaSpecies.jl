
struct ReactionFormula
  subs
  prods
  substoich
  prodstoich
  reverse::Bool
  function ReactionFormula(subs, prods, substoich, prodstoich, reverse)
    subs, substoich = sort_by_species(combine_same(subs, substoich)...)
    prods, prodstoich = sort_by_species(combine_same(prods, prodstoich)...)
    new(subs, prods, substoich, prodstoich, reverse)
  end
end


# Fractional coefficients ("0.5O2[X]") are common in wall/surface reactions,
# so a decimal prefix parses as Float64; a plain integer prefix stays Int.
function get_stoich(v::AbstractString)
  m = match(r"^(\d+(?:\.\d+)?)", v)
  isnothing(m) && return 1
  return occursin('.', m[1]) ? parse(Float64, m[1]) : parse(Int, m[1])
end
remove_stoich(s::AbstractString) = replace(s, r"^\d+(?:\.\d+)?" => s"")

# Keep stoichiometry vectors concretely typed: all-Int stays Int, any
# fractional coefficient promotes the whole vector to Float64.
normalize_stoich(v) = all(x -> x isa Integer, v) ? Int.(v) : Float64.(v)


function combine_same(sp, stoich)
  sp = copy(sp)
  stoich = copy(stoich)
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

function _show_side(io::IO, species, stoich)
  if !isone(stoich[1])
    show(io, stoich[1])
  end
  show(io, species[1])
  for i in Iterators.drop(eachindex(species), 1)
    print(io, '+')
    if !isone(stoich[i])
      show(io, stoich[i])
    end
    show(io, species[i])
  end
end

function Base.show(io::IO, recipe::ReactionFormula)
  _show_side(io, recipe.subs, recipe.substoich)
  print(io, recipe.reverse ? "<-->" : "-->")
  _show_side(io, recipe.prods, recipe.prodstoich)
end

function sort_by_species(sp, stoich)
  p = sortperm(sp, rev=true)
  return sp[p], stoich[p]
end

function parse_reaction(str)
  str = replace(str, r"\s" => "")
  dir = match(r"(-->|<-->|<--)", str)[1]
  reverse_reaction = dir == "<-->"
  # Split on the arrow that actually matched, so "<-->" does not leave a stray
  # '<' on the left-hand side.
  lhs, rhs = dir == "<--" ? split(str, dir)[end:-1:1] : split(str, dir)
  substrates_str = split(lhs, r"(?<!\[|,)\+(?!,|\]|\s\]|\s,)")
  products_str = split(rhs, r"(?<!\[|,)\+(?!,|\]|\s\]|\s,)")
  substoich = normalize_stoich(get_stoich.(substrates_str))
  prodstoich = normalize_stoich(get_stoich.(products_str))
  subs = Species.(remove_stoich.(substrates_str))
  prods = Species.(remove_stoich.(products_str))
  subs, substoich = sort_by_species(combine_same(subs, substoich)...)
  prods, prodstoich = sort_by_species(combine_same(prods, prodstoich)...)
  return ReactionFormula(
    subs,
    prods,
    substoich,
    prodstoich,
    reverse_reaction
  )
end
ReactionFormula(str) = parse_reaction(str)

macro p_str(s)
  if contains(s, r"(-->|<-->|<--)")
    return parse_reaction(s)
  else
    return Species(s)
  end

end

"""
    ```julia 
    apply_tree(t::SpeciesTree, reactions::ReactionFormula)
    ```

Apply the species tree to a reaction and return an array with tuples containing 
the scaling factor and the new reaction. 
This means for each participating species it is checked, if it is a leaf of the tree. 
If not, a new reaction is created for each descendent species. 
If the non-leaf species is a product, the branching factor is one over the number of 
descendants (effectively assuming that the total reaction rate is distributed equally
over all possible products).

"""
function apply_tree(t::SpeciesTree, reaction::ReactionFormula)
  sub_leaves = [matching_leaves(t, s) for s in reaction.subs]
  prod_leaves = [matching_leaves(t, s) for s in reaction.prods]

  # A species with no match expands to nothing, which would drop the reaction
  # from the chemistry without a trace.
  missing_species = [s for (s, l) in zip([reaction.subs; reaction.prods],
                                         [sub_leaves; prod_leaves]) if isempty(l)]
  if !isempty(missing_species)
    @warn "Dropping reaction: species not found in tree" reaction missing_species
    return Tuple{Float64,ReactionFormula}[]
  end

  substrate_vectors = Iterators.product(sub_leaves...) .|> collect
  products_vectors = Iterators.product(prod_leaves...) .|> collect
  branching_factor = 1 / length(products_vectors)
  [(branching_factor, ReactionFormula(subs, prods, reaction.substoich, reaction.prodstoich, reaction.reverse)) for subs in substrate_vectors for prods in products_vectors]
end

function apply_tree(t::SpeciesTree, reactions::ReactionFormula...)
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
apply_tree(t::SpeciesTree, reactions::Vector{ReactionFormula}) = PlasmaSpecies.apply_tree(t, reactions...)


ismassbalanced(r::ReactionFormula) = sum(mass.(r.prods) .* r.prodstoich) ≈ sum(mass.(r.subs) .* r.substoich)
ischargebalanced(r::ReactionFormula) = sum(to_value.(charge.(r.prods)) .* r.prodstoich) == sum(to_value.(charge.(r.subs)) .* r.substoich)
isreactive(r::ReactionFormula) = !(r.subs == r.prods && r.substoich == r.prodstoich)


struct PlasmaReaction{R,F<:ReactionFormula}
    rate::R
    formula::F
end

isreactive(r::PlasmaReaction) = isreactive(r.formula)

"""
    scale_rate(rate, factor)

Scale a reaction rate by a branching factor. A number scales directly; anything
else is treated as callable and wrapped, so rate coefficients that depend on
temperature or reduced field survive [`apply_tree`](@ref).
"""
scale_rate(rate::Number, factor) = rate * factor
scale_rate(rate, factor) = (args...) -> rate(args...) * factor

function apply_tree(t::SpeciesTree, pr::PlasmaReaction)
  [PlasmaReaction(scale_rate(pr.rate, branching_factor), formula) for (branching_factor, formula) in apply_tree(t, pr.formula)]
end

function apply_tree(t::SpeciesTree, prs::PlasmaReaction...)
  out = []
  for pr in prs
    try
      append!(out, apply_tree(t, pr))
    catch error
      @warn "Skipping reaction because an error occured while making " pr error
    end
  end
  return out
end
apply_tree(t::SpeciesTree, prs::Vector{<:PlasmaReaction}) = apply_tree(t, prs...)

function thermal_source_term(pr::PlasmaReaction)
  ΔE = reaction_energy(pr)
  isnothing(ΔE) && error("Cannot compute thermal source term: one or more species have unknown energy.")
  stoichs = pr.formula.substoich
  (ns...) -> pr.rate * prod(n^st for (n, st) in zip(ns, stoichs)) * (-ΔE)
end

function reaction_energy(r::PlasmaReaction)
  f = r.formula
  any(isnothing ∘ energy, f.subs) && return nothing
  any(isnothing ∘ energy, f.prods) && return nothing
  sum(energy(sp) * st for (sp, st) in zip(f.prods, f.prodstoich)) -
  sum(energy(sp) * st for (sp, st) in zip(f.subs, f.substoich))
end


