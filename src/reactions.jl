
"""
    ReactionFormula(subs, prods, substoich, prodstoich, reverse)
    ReactionFormula(str::AbstractString)

Stoichiometry of one reaction: substrate and product species with their
coefficients, plus whether it is reversible. No rate — pair it with one using
[`PlasmaReaction`](@ref).

## Fields

- `subs`, `prods`: `Vector{Species}`.
- `substoich`, `prodstoich`: coefficients, `Int` throughout unless a fractional
  one appears, which promotes the whole side to `Float64`.
- `reverse::Bool`: `true` for `<-->`.

The constructor normalises both sides: repeated species are combined
(`e + e` becomes `2e`) and species are sorted, so two formulas written
differently but meaning the same thing compare and print identically.

`==` and `hash` are field-wise over the normalised form (`reverse` included),
so formulas work as `Dict` keys and in `Set`/`unique`/`in`. Direction is *not*
canonicalised: a reversible `A <--> B` is not equal to `B <--> A`.

From a string, see [`parse_reaction`](@ref); the [`@p_str`](@ref) macro is the
short form.
"""
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

# Comparison is field-wise on the *normalised* form, which is why it is this
# simple: the inner constructor has already combined repeated species and
# sorted both sides, so there is no canonicalisation left to do here. A
# reversible `A <--> B` is deliberately *not* equal to `B <--> A` -- the
# constructor does not canonicalise direction, so comparison does not invent
# it either.
Base.:(==)(a::ReactionFormula, b::ReactionFormula) =
  a.subs == b.subs &&
  a.prods == b.prods &&
  a.substoich == b.substoich &&
  a.prodstoich == b.prodstoich &&
  a.reverse == b.reverse

# Must mirror `==` exactly (same field set), for the same reason as
# `hash(::Species)`: Julia's `Dict`/`Set` require `a == b` to imply
# `hash(a) == hash(b)`. The fields here are `Vector`s -- fresh objects on every
# construction -- so the default `objectid`-derived hash makes two separately
# built formulas that mean the same thing hash differently, silently breaking
# `Dict{ReactionFormula,V}`, `Set`, `unique` and `in`. `0x9c4b2f7a...` is just a
# fixed salt distinguishing this hash from an unrelated type that happened to
# hash its own fields the same way; any constant works.
const _REACTIONFORMULA_HASH_SALT = UInt === UInt64 ? 0x9c4b2f7ae1d63508 : 0x9c4b2f7a
function Base.hash(r::ReactionFormula, h::UInt)
  h = hash(r.subs, h)
  h = hash(r.prods, h)
  h = hash(r.substoich, h)
  h = hash(r.prodstoich, h)
  h = hash(r.reverse, h)
  return hash(_REACTIONFORMULA_HASH_SALT, h)
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

"""
    parse_reaction(str) -> ReactionFormula

Parse a reaction string in LoKI-B notation. Sides are separated by an arrow and
terms by `+`:

    e + N2 --> 2e + N2[+]
    N2[X,vib=0] <--> N2[X,vib=1]     # reversible
    N2[X,vib=1] <-- N2[X,vib=0]      # normalised to the forward direction

Each term is an optional stoichiometric coefficient followed by a species
string (see [`Species(::String)`](@ref)). A decimal coefficient
(`0.5O2[X]`, common in wall reactions) makes that side `Float64`. Whitespace is
ignored, and a `+` inside brackets (`N2[+]`) is not a term separator.
"""
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

"""
    p"..."

Parse a species or a reaction, whichever the string looks like: it becomes a
[`ReactionFormula`](@ref) if it contains an arrow, otherwise a
[`Species`](@ref).

```julia
p"N2[X,vib=1]"            # Species
p"e + N2 --> 2e + N2[+]"  # ReactionFormula
```
"""
macro p_str(s)
  if contains(s, r"(-->|<-->|<--)")
    return parse_reaction(s)
  else
    return Species(s)
  end

end

"""
    apply_tree(t::SpeciesTree, r::ReactionFormula)  -> Vector{Tuple{Float64,ReactionFormula}}
    apply_tree(t::SpeciesTree, pr::PlasmaReaction)  -> Vector{PlasmaReaction}

Expand a reaction written over lumped species into one reaction per combination
of actual states in `t`.

Every species in the formula is resolved to the leaves of `t` it matches, via
[`matching_leaves`](@ref) — an exact state matches itself, a lumped one
(`N2[X]`) matches its whole subtree, and a range (`vib=0-4`) matches the levels
it contains. The result is one reaction per element of the product of those
sets.

Each gets a branching factor of `1 / (number of product combinations)`: the
total rate is split equally over the possible outcomes. For a
[`ReactionFormula`](@ref) the factor is returned alongside the reaction; for a
[`PlasmaReaction`](@ref) it is folded into the rate by [`scale_rate`](@ref).

Reactions whose species are not in the tree at all are dropped with a warning,
rather than silently expanding to nothing. Passing several reactions (varargs or
a vector) concatenates the results and skips any that throw, again with a
warning.

```julia
julia> t = SpeciesTree(["e", "N2[X,vib=0]", "N2[X,vib=1]", "N2[X,vib=2]"]);

julia> apply_tree(t, p"e + N2[X,vib=0-1] --> e + N2[X,vib=2]")
2-element Vector{Tuple{Float64, ReactionFormula}}:
 (1.0, e+N2[X,vib=0]-->e+N2[X,vib=2])
 (1.0, e+N2[X,vib=1]-->e+N2[X,vib=2])
```
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


"""
    ismassbalanced(r::ReactionFormula) -> Bool

Whether both sides have the same total mass, to within `≈`. Requires every
species to have a resolvable mass — errors on a [`StringGas`](@ref). Not
exported.
"""
ismassbalanced(r::ReactionFormula) = sum(mass.(r.prods) .* r.prodstoich) ≈ sum(mass.(r.subs) .* r.substoich)

"""
    ischargebalanced(r::ReactionFormula) -> Bool

Whether both sides carry the same total charge. Not exported.
"""
ischargebalanced(r::ReactionFormula) = sum(to_value.(charge.(r.prods)) .* r.prodstoich) == sum(to_value.(charge.(r.subs)) .* r.substoich)

"""
    isreactive(r) -> Bool

Whether a [`ReactionFormula`](@ref) or [`PlasmaReaction`](@ref) changes
anything: `false` when both sides hold the same species with the same
coefficients. [`apply_tree`](@ref) can produce such identities when a lumped
reaction expands onto a single state, and they contribute nothing to a
chemistry.
"""
isreactive(r::ReactionFormula) = !(r.subs == r.prods && r.substoich == r.prodstoich)


"""
    PlasmaReaction(rate, formula)

A [`ReactionFormula`](@ref) with a rate attached.

`rate` is deliberately untyped: a number is a constant rate coefficient,
anything callable is evaluated by the consumer (a temperature- or
field-dependent fit, an interpolation from a Boltzmann solver). Only
[`scale_rate`](@ref) touches it here.

For that reason `PlasmaReaction` has no `==`/`hash` of its own and falls back
to identity: a `rate` is often a closure, and closures compare by identity, so
equality of a whole reaction is not well-defined. Compare the `formula` fields
instead -- [`ReactionFormula`](@ref) has proper `==`/`hash`.
"""
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

"""
    thermal_source_term(pr::PlasmaReaction) -> Function

Build `(n₁, n₂, ...) -> rate * ∏ nᵢ^stoichᵢ * (-ΔE)`, the power released into
the gas by `pr` at the given substrate densities. Arguments are the substrate
densities in the order of `pr.formula.subs`.

Sign convention: an exothermic reaction (`ΔE < 0`, see
[`reaction_energy`](@ref)) gives a positive source. Errors unless every species
in the reaction has an `energy`; set them with [`apply_energy!`](@ref).

Assumes a constant `rate` — the product is formed directly, without evaluating
a callable rate.
"""
function thermal_source_term(pr::PlasmaReaction)
  ΔE = reaction_energy(pr)
  isnothing(ΔE) && error("Cannot compute thermal source term: one or more species have unknown energy.")
  stoichs = pr.formula.substoich
  (ns...) -> pr.rate * prod(n^st for (n, st) in zip(ns, stoichs)) * (-ΔE)
end

"""
    reaction_energy(r::PlasmaReaction) -> Union{Number,Nothing}

Energy balance `Σ E(products) - Σ E(substrates)`, weighted by stoichiometry.
Negative for an exothermic reaction. Returns `nothing` if any species has no
`energy` set, which is the expected state for a chemistry that has not had
[`apply_energy!`](@ref) applied.

The unit is whatever [`energy`](@ref) carries.
"""
function reaction_energy(r::PlasmaReaction)
  f = r.formula
  any(isnothing ∘ energy, f.subs) && return nothing
  any(isnothing ∘ energy, f.prods) && return nothing
  sum(energy(sp) * st for (sp, st) in zip(f.prods, f.prodstoich)) -
  sum(energy(sp) * st for (sp, st) in zip(f.subs, f.substoich))
end


