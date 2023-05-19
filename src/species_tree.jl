using AbstractTrees

struct SpeciesTree
    x::Union{Species,Nothing}
    children::Vector{SpeciesTree}
end
export SpeciesTree

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


insert_child!(t::SpeciesTree, child::SpeciesTree) = push!(children(t), child)
export insert_child!

#TODO:finish implementation, make tests pass
function combine(t1::SpeciesTree, t2::SpeciesTree)
    return SpeciesTree(nodevalue(t1), [children(t1);children(t2)])
end

