struct Species
    gas::String
    charge::Union{String,Nothing}
    electronic_state::Union{String,Nothing}
    vibrational_state::Union{String,Nothing}
    rotational_state::Union{String,Nothing}
end

function Species(str::String)
    regex = r"(^(?:\w|\d)*)(?:\(([+\-0]*)?,?([^,\)]*)?(?:,v=([^,\)]*)(?:,J=([^,\)]*))?)?\))?"
    m = match(regex, str)
    Species(m...)
end


gas(sp::Species) = sp.gas
charge(sp::Species) = sp.charge
electronic_state(sp::Species) = sp.electronic_state
vibrational_state(sp::Species) = sp.vibrational_state
rotational_state(sp::Species) = sp.rotational_state
export Species, gas, electronic_state, vibrational_state, rotational_state

function get_parent_species(sp::Species)
    if !isnothing(rotational_state(sp))
        return Species(
            gas(sp),
            charge(sp),
            electronic_state(sp),
            vibrational_state(sp),
            nothing
        )
    elseif !isnothing(vibrational_state(sp))
        return Species(
            gas(sp),
            charge(sp),
            electronic_state(sp),
            nothing,
            nothing
        )
    elseif !isnothing(electronic_state(sp))
        return Species(
            gas(sp),
            charge(sp),
            nothing,
            nothing,
            nothing
        )
    elseif !isnothing(charge(sp))
        return Species(
            gas(sp),
            nothing,
            nothing,
            nothing,
            nothing
        )
    else
        return nothing
    end
end
export get_parent_species

function is_parent_species(sp::Species, parent::Species)
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

Base.show(io::IO, sp::Species) =
    Base.show(io,
        gas(sp) *
        (isnothing(electronic_state(sp)) ? "" : "(" * electronic_state(sp) *
         (isnothing(vibrational_state(sp)) ? ")" : ",v=" * vibrational_state(sp) *
          (isnothing(rotational_state(sp)) ? ")" : ",J=" * rotational_state(sp) * ")")
         )
        )
    )
Base.print(io::IO, sp::Species) =
    Base.print(io,
        gas(sp) *
        (isnothing(electronic_state(sp)) ? "" : "(" * electronic_state(sp) *
         (isnothing(vibrational_state(sp)) ? ")" : ",v=" * vibrational_state(sp) *
          (isnothing(rotational_state(sp)) ? ")" : ",J=" * rotational_state(sp) * ")")
         )
        )
    )

