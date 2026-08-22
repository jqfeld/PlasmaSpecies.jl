```@meta
CurrentModule = PlasmaSpecies
```

# Species trees

A [`SpeciesTree`](@ref) holds the parent/child hierarchy of species: each child
is a more resolved state of its parent, obtained by filling in one more field.
Bulk `N2` has the child `N2[X]`, which has the children `N2[X,vib=0]`,
`N2[X,vib=1]`, and so on.

The tree exists so that a chemistry can be *written* at one level of detail and
*solved* at another. Reactions are written over whatever species the data
source names — often lumped, `e + N2[X] --> ...` — while the model integrates
the leaves. [`apply_tree`](@ref) bridges the two.

## Building

Constructing from a species inserts its whole ancestor chain, found by repeated
[`get_parent_species`](@ref). Constructing from several merges those chains,
sharing every common ancestor. The root is a node holding `nothing`.

```@repl trees
using PlasmaSpecies
t = SpeciesTree(["e", "N2[X,vib=0]", "N2[X,vib=1]", "N2[X,vib=2]", "N2[+]"])
```

Accepted arguments: a `Species` or `String`, a vector of either, or varargs.
Add to an existing tree with `push!`, or combine two with `merge!`:

```@repl trees
push!(t, "N2[X,vib=3]")
```

`t[sp] = subtree` grafts a separately built branch onto an existing tree,
replacing the node for `sp` value and children alike.

## Querying

```@repl trees
leaves(t)
t["N2[X]"]
matching_leaves(t, Species("N2[X,vib=0-1]"))
matching_leaves(t, Species("O2"))
```

[`leaves`](@ref) gives the most resolved states — the model's unknowns.
Indexing finds one node and returns its subtree, stopping at the first match.
[`matching_leaves`](@ref) is the pattern query: it descends every branch and
collects the leaves under each node that the given species matches, so a
wildcard field or a `vib=0-4` range expands to every level actually present.

The tree implements `AbstractTrees`, so `print_tree`, `Leaves` and the rest of
that interface work directly on it.

## Energies

[`apply_energy!`](@ref) fills the `energy` field of every leaf from a function
of its state labels:

```@repl trees
tv = SpeciesTree(["N2[X,vib=0]", "N2[X,vib=1]", "N2[X,vib=2]"])
apply_energy!(tv, q -> 0.29 * q.vib)   # eV, harmonic
energy.(leaves(tv))
```

The callback receives a named tuple `(elec, vib, rot)`; the unit is the
caller's choice, and only [`reaction_energy`](@ref) and
[`thermal_source_term`](@ref) consume it. Leaves are rebuilt with
[`with_fields`](@ref), so every other field survives.

```@docs; canonical=false
SpeciesTree
leaves
matching_leaves
apply_energy!
```

`push!`, `merge!`, `getindex` and `setindex!` are extended for `SpeciesTree` as
described above; their docstrings are available at the REPL.
