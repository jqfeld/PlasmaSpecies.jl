using AbstractTrees

using Base: setindex
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

Base.push!(t::SpeciesTree, x) = Base.merge!(t, SpeciesTree(x))

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

