using AbstractTrees

using Base: setindex
"""
    SpeciesTree(sp::Species)
    SpeciesTree(str::String)
    SpeciesTree(species)          # Vector or varargs of Species / String

The species hierarchy: a node holds one [`Species`](@ref) (or `nothing` for the
root) and its children, each a more resolved state of its parent.

Building from a species inserts its whole ancestor chain, obtained by repeated
[`get_parent_species`](@ref), so `"N2[X,vib=1]"` alone yields
`nothing → N2 → N2[X] → N2[X,vib=1]`. Building from several species merges those
chains, sharing every common ancestor.

```julia
julia> SpeciesTree(["e", "N2[X,vib=0]", "N2[X,vib=1]", "N2[+]"])
nothing
├─ e
├─ N2
│  └─ N2[X]
│     ├─ N2[X,vib=0]
│     └─ N2[X,vib=1]
└─ N2[+]
```

The tree implements `AbstractTrees`, so `print_tree`, `Leaves` and friends work
on it. Supported operations: `push!` and `merge!` to add species, `t[sp]` /
`t["N2[X]"]` to fetch a subtree (`nothing` if absent), `t[sp] = subtree` to
replace one, [`leaves`](@ref) and [`matching_leaves`](@ref) to collect states,
and [`apply_energy!`](@ref) to fill in level energies.

The leaves are the states a model actually integrates; interior nodes are the
lumped species reactions may be written over. [`apply_tree`](@ref) uses that to
expand a reaction written for `N2[X]` into one per vibrational level.
"""
mutable struct SpeciesTree
    x::Union{Species,Nothing}
    children::Vector{SpeciesTree}
end

Base.show(io::IO, t::SpeciesTree) = AbstractTrees.print_tree(io,t)
AbstractTrees.children(t::SpeciesTree) = t.children
AbstractTrees.nodevalue(t::SpeciesTree) = t.x

function SpeciesTree(sp::Species)
    parents = []
    t = SpeciesTree(sp, [])
    tmp = sp
    while !isnothing(get_parent_species(tmp))
        push!(parents, get_parent_species(tmp))
        tmp = get_parent_species(tmp)
    end
    while !isempty(parents)
        t = SpeciesTree(
            popfirst!(parents),
            [t]
        )
    end
    return SpeciesTree(nothing, [t])
end
SpeciesTree(str::String) = SpeciesTree(Species(str))
SpeciesTree(str::Vector{String}) = merge!(SpeciesTree.(str)...)
SpeciesTree(str::String...) = SpeciesTree(String[str...])
SpeciesTree(sps::Vector{<:Species}) = merge!(SpeciesTree.(sps)...)
SpeciesTree(sps::Species...) = SpeciesTree(Species[sps...])

"""
    merge!(left::SpeciesTree, right::SpeciesTree...) -> SpeciesTree

Merge `right` into `left` in place, recursing wherever both hold the same
species so shared ancestors are not duplicated. Returns `left`.
"""
function Base.merge!(left::SpeciesTree, right::SpeciesTree)
    for nv in children(right)
        if nodevalue(nv) ∈ nodevalue.(children(left))
            merge!(children(left)[findfirst(==(nodevalue(nv)), nodevalue.(children(left)))], nv)
        else
            push!(children(left), nv)
        end
    end
    return left
end

function Base.merge!(left::SpeciesTree, others::SpeciesTree...)
    reduce(merge!, [left; others...])
end

"""
    push!(t::SpeciesTree, x) -> SpeciesTree

Add a species (or string) and its ancestor chain to `t`. Equivalent to
`merge!(t, SpeciesTree(x))`.
"""
Base.push!(t::SpeciesTree, x) = Base.merge!(t, SpeciesTree(x))

"""
    t[sp::Species] -> Union{SpeciesTree,Nothing}
    t[str::String] -> Union{SpeciesTree,Nothing}

Depth-first search for the node holding `sp`, returned as a subtree, or
`nothing` if it is not in `t`. Stops at the first match — use
[`matching_leaves`](@ref) to collect every state a wildcard or range matches.
"""
function Base.getindex(t::SpeciesTree, s::Species)
    for c in children(t)
        if nodevalue(c) == s
            return c
        end

        # if not found on this level, recursively ascend the tree
        out = Base.getindex(c, s)
        if !isnothing(out)
            return out
        end
    end
end

Base.getindex(t::SpeciesTree, s::String) = t[Species(s)]

"""
    t[sp] = subtree

Replace the node holding `sp` in `t` with the node holding `sp` in `subtree`,
value and children alike. Throws `KeyError` if `sp` is missing from either.
Used to graft a separately built branch (e.g. a full vibrational ladder) onto
an existing tree.
"""
function Base.setindex!(t::SpeciesTree, tn::SpeciesTree, s::Species)
    target = t[s]
    isnothing(target) && throw(KeyError(s))
    source = tn[s]
    isnothing(source) && throw(KeyError(s))
    target.x = source.x
    target.children = source.children
    return t
end

Base.setindex!(t::SpeciesTree, tn::SpeciesTree, s::String) = setindex!(t, tn, Species(s))

"""
    leaves(t::SpeciesTree) -> Vector{Species}

Every leaf species of `t`, in tree order. These are the most resolved states in
the hierarchy — the ones a model carries as unknowns.
"""
leaves(t::SpeciesTree) = Leaves(t) |> collect .|> nodevalue

"""
    matching_leaves(t::SpeciesTree, formula::Species) -> Vector{Species}

Collect every leaf of `t` matched by `formula` (via `nodevalue(c) in formula`,
see [`Base.in`](@ref)/[`species_matches`](@ref)), descending into every branch
rather than stopping at the first match like `t[formula]`/`Base.getindex`
does. Whenever `formula` matches a node — exactly, or via a `UnitRange`
vibrational/rotational field containing that node's level — the whole
subtree's leaves are taken and that branch is not searched further; otherwise
the search continues into the node's children. This is what lets a reaction
formula species written with a vibrational range (e.g. `vib=0-4`) expand, in
[`apply_tree`](@ref), to every matching level actually present in the tree,
even when those levels are separate sibling nodes.
"""
function matching_leaves(t::SpeciesTree, formula::Species)
  result = Species[]
  for c in children(t)
    if nodevalue(c) in formula
      append!(result, leaves(c))
    else
      append!(result, matching_leaves(c, formula))
    end
  end
  return result
end



"""
    apply_energy!(t::SpeciesTree, energy) -> SpeciesTree

Set the `energy` field of every leaf of `t` in place. `energy` is called with a
named tuple `(elec, vib, rot)` of that leaf's state labels and returns its level
energy; the unit is the caller's choice.

```julia
apply_energy!(t, q -> ωe * (q.vib + 0.5) - ωeχe * (q.vib + 0.5)^2)
```

Leaves are replaced via [`with_fields`](@ref), so every other field survives.
Interior nodes are left alone. [`reaction_energy`](@ref) consumes what this
sets.
"""
function apply_energy!(tr, energy)
    for l in AbstractTrees.Leaves(tr)
        species = l.x
        l.x = with_fields(
            species,
            energy = energy(
                (
                    elec = electronic_state(species),
                    vib = vibrational_state(species),
                    rot = rotational_state(species),
                )
            )
        )
    end
    return tr
end

