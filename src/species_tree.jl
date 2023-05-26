using AbstractTrees

struct SpeciesTree
    x::Union{SpeciesNode,Nothing}
    children::Vector{SpeciesTree}
end
export SpeciesTree

AbstractTrees.children(t::SpeciesTree) = t.children
AbstractTrees.nodevalue(t::SpeciesTree) = t.x

function SpeciesTree(sp::SpeciesNode)
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
SpeciesTree(str::String) = SpeciesTree(SpeciesNode(str))
SpeciesTree(str::Vector{String}) = merge!(SpeciesTree.(str)...)

insert_child!(t::SpeciesTree, child::SpeciesTree) = push!(children(t), child)
export insert_child!

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



function Base.getindex(t::SpeciesTree, s::String)
    node = SpeciesNode(s)
    for c in children(t)
        if nodevalue(c) == node
            return c
        end
        out = Base.getindex(c, s)
        if !isnothing(out)
            return out
        end
    end
end
