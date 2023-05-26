struct SpeciesNode
    gas::String
    charge::Union{String,Nothing}
    electronic_state::Union{String,Nothing}
    vibrational_state::Union{String,Nothing}
    rotational_state::Union{String,Nothing}
end

function SpeciesNode(str::String)
    regex = r"(^(?:\w|\d)*)(?:\(([+\-]+)?,?([^,\)]+)?(?:,v=([^,\)]*)(?:,J=([^,\)]*))?)?\))?"
    m = match(regex, str)
    SpeciesNode(m...)
end


gas(sp::SpeciesNode) = sp.gas
charge(sp::SpeciesNode) = sp.charge
electronic_state(sp::SpeciesNode) = sp.electronic_state
vibrational_state(sp::SpeciesNode) = sp.vibrational_state
rotational_state(sp::SpeciesNode) = sp.rotational_state

export SpeciesNode, gas, electronic_state, vibrational_state, rotational_state

function get_parent_species(sp::SpeciesNode)
    if !isnothing(rotational_state(sp))
        return SpeciesNode(
            gas(sp),
            charge(sp),
            electronic_state(sp),
            vibrational_state(sp),
            nothing
        )
    elseif !isnothing(vibrational_state(sp))
        return SpeciesNode(
            gas(sp),
            charge(sp),
            electronic_state(sp),
            nothing,
            nothing
        )
    elseif !isnothing(electronic_state(sp))
        return SpeciesNode(
            gas(sp),
            charge(sp),
            nothing,
            nothing,
            nothing
        )
    # elseif !isnothing(charge(sp))
    #     return SpeciesNode(
    #         gas(sp),
    #         nothing,
    #         nothing,
    #         nothing,
    #         nothing
    #     )
    else
        return nothing
    end
end
export get_parent_species

function is_parent_species(sp::SpeciesNode, parent::SpeciesNode)
    if parent != sp && gas(parent) == gas(sp)
        if isnothing(electronic_state(parent)) || electronic_state(parent) == electronic_state(sp)
            return isnothing(vibrational_state(parent)) || vibrational_state(parent) == vibrational_state(sp) ? true : false
        else
            false
        end
    else
        return false
    end
end
export is_parent_species


function Base.print(io::IO, sp::SpeciesNode)
    out = gas(sp)
    open_bracket = false
    if !isnothing(charge(sp))
        open_bracket = true
        out *= "($(charge(sp))"
    end
    if !isnothing(electronic_state(sp))
        out *= !open_bracket ? "(" : ","
        open_bracket = true
        out *= electronic_state(sp)
    end
    if isnothing(vibrational_state(sp))
        out *= open_bracket ? ")" : ""
        return Base.print(io, out)
    else
        out *= ",v=" * vibrational_state(sp)
    end
    if isnothing(rotational_state(sp))
        out *= open_bracket ? ")" : ""
        return Base.print(io, out)
    else
        out *= ",J=" * rotational_state(sp) * ")"
    end
    Base.print(io, out)
end

Base.show(io::IO, sp::SpeciesNode) = Base.show(io, string(sp))
