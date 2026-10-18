
"""
    ReactionFormula(subs, prods, substoich, prodstoich, reverse)
    ReactionFormula(str::AbstractString)

Stoichiometry of one reaction: substrate and product species with their
coefficients, plus whether it is reversible. No rate — pair it with one using
[`PlasmaReaction`](@ref).

## Fields

- `subs`, `prods`: `Vector{Species}`.
- `substoich::Vector{S}`, `prodstoich::Vector{P}`: coefficients. Each side is
  typed on its own: `Int` when all its coefficients are integers, otherwise
  `Float64`, so `O[3P] --> 0.5O2[X]` is a `ReactionFormula{Int,Float64}`.
- `reverse::Bool`: `true` for `<-->`.

The constructor normalises both sides: repeated species are combined
(`e + e` becomes `2e`) and species are sorted, so two formulas written
differently but meaning the same thing compare and print identically. It
throws `ArgumentError` if a side has a different number of species and
coefficients.

`==` and `hash` are field-wise over the normalised form (`reverse` included),
so formulas work as `Dict` keys and in `Set`/`unique`/`in`. Direction is *not*
canonicalised: a reversible `A <--> B` is not equal to `B <--> A`.

From a string, see [`parse_reaction`](@ref); the [`@p_str`](@ref) macro is the
short form.
"""
struct ReactionFormula{S<:Real,P<:Real}
  subs::Vector{Species}
  prods::Vector{Species}
  substoich::Vector{S}
  prodstoich::Vector{P}
  reverse::Bool
  function ReactionFormula(subs, prods, substoich, prodstoich, reverse)
    length(subs) == length(substoich) && length(prods) == length(prodstoich) ||
      throw(ArgumentError("a reaction side needs one coefficient per species"))
    S = _stoich_type(substoich)
    P = _stoich_type(prodstoich)
    subs, substoich = sort_by_species(combine_same(Species[subs...], S[substoich...])...)
    prods, prodstoich = sort_by_species(combine_same(Species[prods...], P[prodstoich...])...)
    new{S,P}(subs, prods, substoich, prodstoich, reverse)
  end
end

# A side's coefficients are Int when they all are, else Float64: concrete, and
# an integer side stays exact next to a fractional one.
_stoich_type(v) = all(x -> x isa Integer, v) ? Int : Float64

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

# Split one side of a reaction into terms at every '+' outside brackets, so any
# charge ("N2[+]", "Ar[+++]", "N2[+,B]") stays inside its species.
function split_terms(side::AbstractString)
  terms = SubString{String}[]
  s = String(side)
  depth = 0
  start = firstindex(s)
  for i in eachindex(s)
    c = s[i]
    if c == '['
      depth += 1
    elseif c == ']'
      depth -= 1
    elseif c == '+' && depth == 0
      push!(terms, SubString(s, start, prevind(s, i)))
      start = nextind(s, i)
    end
  end
  push!(terms, SubString(s, start, lastindex(s)))
  return terms
end

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

# A coefficient of 1 is implicit, and an integral Float64 prints as an integer
# so a side written "2.0N2 + 0.5O2" reads back as "2N2+0.5O2".
function _show_coefficient(io::IO, st)
  isone(st) && return
  isinteger(st) ? print(io, Int(st)) : show(io, st)
end

function _show_side(io::IO, species, stoich)
  _show_coefficient(io, stoich[1])
  show(io, species[1])
  for i in Iterators.drop(eachindex(species), 1)
    print(io, '+')
    _show_coefficient(io, stoich[i])
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
string (see [`Species(::AbstractString)`](@ref)). A decimal coefficient
(`0.5O2[X]`, common in wall reactions) makes that side's coefficients `Float64`;
the other side keeps `Int` (see [`ReactionFormula`](@ref)). Whitespace is
ignored, and a `+` inside brackets (`N2[+]`) is not a term separator.

Leading digits are always the coefficient. A mass number is written in braces
in front of its symbol (see [`parse_formula`](@ref)), so `2{16}O{1}H` is two
¹⁶O¹H and `0.5{14}N2` is half a ¹⁴N₂.
"""
function parse_reaction(str)
  str = replace(str, r"\s" => "")
  arrows = collect(eachmatch(r"<-->|-->|<--", str))
  length(arrows) == 1 || throw(ArgumentError(
    length(arrows) == 0 ?
      "no reaction arrow in $(repr(str)); expected one of -->, <--, <-->" :
      "more than one reaction arrow in $(repr(str))"))
  dir = arrows[1].match
  reverse_reaction = dir == "<-->"
  # Split on the arrow that actually matched, so "<-->" does not leave a stray
  # '<' on the left-hand side.
  lhs, rhs = dir == "<--" ? split(str, dir)[end:-1:1] : split(str, dir)
  substrates_str = split_terms(lhs)
  products_str = split_terms(rhs)
  # The constructor combines repeated species, sorts, and types each side.
  return ReactionFormula(
    Species.(remove_stoich.(substrates_str)),
    Species.(remove_stoich.(products_str)),
    get_stoich.(substrates_str),
    get_stoich.(products_str),
    reverse_reaction,
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
a vector) concatenates the results. Any other error aborts the whole batch: it
is logged with the reaction that caused it and then rethrown, so a chemistry is
never silently built with a reaction missing.

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
  Tuple{Float64,ReactionFormula}[
    (branching_factor, ReactionFormula(subs, prods, reaction.substoich, reaction.prodstoich, reaction.reverse))
    for subs in substrate_vectors for prods in products_vectors]
end

# Run `f()` for one reaction of a batch. An error is logged with the reaction
# it came from and rethrown unchanged, so the batch aborts instead of quietly
# losing a reaction.
function _expanding(f, reaction)
  try
    f()
  catch
    @error "Error while expanding reaction" reaction
    rethrow()
  end
end

function apply_tree(t::SpeciesTree, reactions::ReactionFormula...)
  out = Tuple{Float64,ReactionFormula}[]
  for r in reactions
    append!(out, _expanding(() -> apply_tree(t, r), r))
  end
  return out
end
apply_tree(t::SpeciesTree, reactions::AbstractVector{<:ReactionFormula}) = apply_tree(t, reactions...)


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
    PlasmaReaction(rate, formula; partition = nothing)

A [`ReactionFormula`](@ref) with a rate attached, and optionally how its released energy is shared:
`partition = (; electrons, gas, radiation)`, shares summing to 1, for what the formula does not show (an excimer's
decay `He2[*] --> 2He[X]` radiates: `(radiation = 1,)`; omitted shares are 0). `nothing`: the rule of
[`energy_partition`](@ref). Shares, not energies: [`apply_tree`](@ref) splits a reaction into states of other
energies and passes the shares on.

`rate` is deliberately untyped: a number is a constant rate coefficient,
anything callable is evaluated by the consumer (a temperature- or
field-dependent fit, an interpolation from a Boltzmann solver). Only
[`scale_rate`](@ref) touches it here.

For that reason `PlasmaReaction` has no `==`/`hash` of its own and falls back
to identity: a `rate` is often a closure, and closures compare by identity, so
equality of a whole reaction is not well-defined. Compare the `formula` fields
instead -- [`ReactionFormula`](@ref) has proper `==`/`hash`.
"""
struct PlasmaReaction{R,F<:ReactionFormula,P}
    rate::R
    formula::F
    partition::P
end
PlasmaReaction(rate, formula::ReactionFormula; partition = nothing) = PlasmaReaction(rate, formula, energy_shares(partition))

"""
    energy_shares(p) -> NamedTuple or nothing

A [`PlasmaReaction`](@ref)'s `partition` checked: `nothing`, or `(; electrons, gas, radiation)` (omitted ones 0), each
in [0, 1], summing to 1.
"""
energy_shares(::Nothing) = nothing
function energy_shares(p)
  keys(p) ⊆ (:electrons, :gas, :radiation) ||
    throw(ArgumentError("partition: the shares are `electrons`, `gas`, `radiation`; got $(join(keys(p), ", "))"))
  sh = (electrons = Float64(get(p, :electrons, 0.0)), gas = Float64(get(p, :gas, 0.0)), radiation = Float64(get(p, :radiation, 0.0)))
  all(x -> 0 <= x <= 1, values(sh)) && isapprox(sum(values(sh)), 1; atol = 1e-12) ||
    throw(ArgumentError("partition: shares in [0, 1] summing to 1, got $sh"))
  return sh
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
  [PlasmaReaction(scale_rate(pr.rate, branching_factor), formula, pr.partition) for (branching_factor, formula) in apply_tree(t, pr.formula)]
end

function apply_tree(t::SpeciesTree, prs::PlasmaReaction...)
  out = PlasmaReaction[]
  for pr in prs
    append!(out, _expanding(() -> apply_tree(t, pr), pr))
  end
  return out
end
apply_tree(t::SpeciesTree, prs::AbstractVector{<:PlasmaReaction}) = apply_tree(t, prs...)

"""
    energy_partition(pr::PlasmaReaction; consumed = 0) -> (; electrons, gas, radiation) or nothing

Where the energy a reaction releases per event goes: `Q = −ΔE` ([`reaction_energy`](@ref); negative: absorbed) plus
the kinetic energy `consumed` [the unit of the levels] of each electron it consumes (recombination, attachment: the
levels do not hold it), split by the reaction's `partition`, else by the rule:

- an electron among the products: the electrons take it (Penning, superelastic, associative and stepwise ionisation;
  a threshold they pay);
- else one heavy product: radiation (a single body cannot take the excess as motion; two-body radiative
  recombination, radiative association);
- else the gas (two or more heavy products share it as motion: VT relaxation, quenching, charge transfer, ion
  conversion, dissociative recombination and attachment, three-body recombination on a heavy third body).

`nothing` when a species' level is unknown.
"""
function energy_partition(pr::PlasmaReaction; consumed = 0)
  ΔE = reaction_energy(pr)
  isnothing(ΔE) && return nothing
  total = -ΔE + consumed_electrons(pr.formula) * consumed
  sh = isnothing(pr.partition) ? rule_shares(pr.formula) : pr.partition
  return (electrons = sh.electrons * total, gas = sh.gas * total, radiation = sh.radiation * total)
end
_is_electron(sp) = gas(sp) isa Electron
"The electrons a reaction consumes: those among its substrates less those among its products (0 if it makes them)."
consumed_electrons(f::ReactionFormula) =
  max(0, sum((st for (sp, st) in zip(f.subs, f.substoich) if _is_electron(sp)); init = 0) -
         sum((st for (sp, st) in zip(f.prods, f.prodstoich) if _is_electron(sp)); init = 0))
"The rule's shares of [`energy_partition`](@ref)."
function rule_shares(f::ReactionFormula)
  any(_is_electron, f.prods) && return (electrons = 1.0, gas = 0.0, radiation = 0.0)
  heavy = sum((st for (sp, st) in zip(f.prods, f.prodstoich) if !_is_electron(sp)); init = 0)
  heavy == 1 && return (electrons = 0.0, gas = 0.0, radiation = 1.0)
  return (electrons = 0.0, gas = 1.0, radiation = 0.0)
end

"""
    energy_source_terms(pr::PlasmaReaction) -> Function

Build `(n₁, n₂, ...; rate_args=(), consumed=0) -> (; electrons, gas, radiation)`, the power the reaction gives each
channel at the given substrate densities: `rate * ∏ nᵢ^stoichᵢ` times the shares of [`energy_partition`](@ref)
(`consumed`: a consumed electron's kinetic energy, the local mean energy of a fluid). Positional arguments are the
substrate densities in the order of `pr.formula.subs`; the function throws `ArgumentError` if their number differs.
A numeric `rate` is used as is and `rate_args` is ignored; a callable `rate` is evaluated as `rate(rate_args...)`.
Errors unless every species in the reaction has an `energy` (an electron's is 0).
"""
function energy_source_terms(pr::PlasmaReaction)
  isnothing(reaction_energy(pr)) && error("Cannot compute the energy source terms: one or more species have unknown energy.")
  stoichs = pr.formula.substoich
  rate = pr.rate
  return function (ns...; rate_args=(), consumed=0)
    length(ns) == length(stoichs) || throw(ArgumentError(
      "expected $(length(stoichs)) substrate densities for $(pr.formula), got $(length(ns))"))
    r = _evaluate_rate(rate, rate_args) * prod(n^st for (n, st) in zip(ns, stoichs))
    q = energy_partition(pr; consumed)
    return (electrons = r * q.electrons, gas = r * q.gas, radiation = r * q.radiation)
  end
end

"""
    thermal_source_term(pr::PlasmaReaction) -> Function

Build `(n₁, n₂, ...; rate_args=(), consumed=0) -> power`, the power released into the gas by `pr` at the given
substrate densities: the `gas` channel of [`energy_source_terms`](@ref), `rate * ∏ nᵢ^stoichᵢ` times the gas's share
of the released energy ([`energy_partition`](@ref)). Positional arguments are the substrate densities in the order of
`pr.formula.subs`; the function throws `ArgumentError` if their number differs.

A numeric `rate` is used as is and `rate_args` is ignored. A callable `rate` is
evaluated as `rate(rate_args...)` on every call, so the closure is built once
and can follow a changing electron temperature:

```julia
f = thermal_source_term(pr)          # pr.rate = Te -> 1e-15 * sqrt(Te)
f(n_e, n_N2; rate_args=(Te,))
```

Rates wrapped by [`apply_tree`](@ref) (see [`scale_rate`](@ref)) take the same
`rate_args` as the original.

Sign convention: an exothermic reaction (`ΔE < 0`, see
[`reaction_energy`](@ref)) gives a positive source where the gas takes it. Errors unless every species
in the reaction has an `energy` (an electron's is 0); set them with [`apply_energy!`](@ref).
"""
function thermal_source_term(pr::PlasmaReaction)
  isnothing(reaction_energy(pr)) && error("Cannot compute thermal source term: one or more species have unknown energy.")
  f = energy_source_terms(pr)
  return (ns...; rate_args=(), consumed=0) -> f(ns...; rate_args, consumed).gas
end

_evaluate_rate(rate::Number, _) = rate
_evaluate_rate(rate, args) = rate(args...)

"""
    reaction_energy(r::PlasmaReaction) -> Union{Number,Nothing}

Energy balance `Σ E(products) - Σ E(substrates)`, weighted by stoichiometry.
Negative for an exothermic reaction. An electron's level is 0 unless set (it has no internal level; its kinetic
energy is the consumer's: [`energy_partition`](@ref)'s `consumed`). Returns `nothing` if any other species has no
`energy` set, which is the expected state for a chemistry that has not had [`apply_energy!`](@ref) applied.

The unit is whatever [`energy`](@ref) carries.
"""
function reaction_energy(r::PlasmaReaction)
  f = r.formula
  any(isnothing ∘ level, f.subs) && return nothing
  any(isnothing ∘ level, f.prods) && return nothing
  sum(level(sp) * st for (sp, st) in zip(f.prods, f.prodstoich)) -
  sum(level(sp) * st for (sp, st) in zip(f.subs, f.substoich))
end
"A species' level for an energy balance: its `energy`, an electron's 0 when unset."
level(sp) = (E = energy(sp); isnothing(E) && _is_electron(sp) ? 0 : E)


