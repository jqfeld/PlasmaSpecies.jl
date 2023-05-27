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
SpeciesTree(str::Vector{String}) = merge!(SpeciesTree.(str)...)
SpeciesTree(str::String...) = SpeciesTree(String[str...])

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
        out = Base.getindex(c, s)
        if !isnothing(out)
            return out
        end
    end
end

Base.getindex(t::SpeciesTree, s::String) = t[Species(s)]

Base.show(io::IO, t::SpeciesTree) = AbstractTrees.print_tree(io,t)
