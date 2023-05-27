
# NOTE: I am not sure the separation of Positive, Negative and Neutral is
# necessarry. Probably, Charge{Val} would be enough, but maybe with the three
# types some nice dispatch magic might be possible later.
abstract type Charge end
struct Positive <: Charge
    value::Int8
end
struct Negative <: Charge
    value::Int8
end
struct Neutral <: Charge end
Charge(::Nothing) = Neutral()

function Charge(s::AbstractString)
    m = match(r"^(\++|\-+|0?)$", s)
    if begin
        m = match(r"^\++$", s)
        !isnothing(m)
    end
        return Positive(Int8(length(m.match)))
    elseif begin
        m = match(r"^\-+$", s)
        !isnothing(m)
    end
        return Negative(Int8(length(m.match)))
    elseif begin
        m = match(r"^0?$", s)
        !isnothing(m)
    end
        return Neutral()
    else
        error("Not able to parse charge")
    end
end

Base.show(io::IO, c::Positive) = Base.print(io::IO, '+'^c.value)
Base.show(io::IO, c::Negative) = Base.print(io::IO, '-'^c.value)
Base.show(io::IO, ::Neutral) = Base.print(io::IO, "")


struct Species{S,C<:Charge,E,V,J}
    gas::S
    charge::C
    electronic_state::E
    vibrational_state::V
    rotational_state::J
end
Base.:(==)(x::S, y::S) where {S<:Species} =
    gas(x) == gas(y) &&
    charge(x) == charge(y) &&
    electronic_state(x) == electronic_state(y) &&
    vibrational_state(x) == vibrational_state(y) &&
    rotational_state(x) == rotational_state(y)

function Species(str::String)
    regex = r"(^(?:\w|\d)*)(?:\(([+\-]+)?,?([^,\)]+)?(?:,v=([^,\)]*)(?:,J=([^,\)]*))?)?\))?"
    m = collect(match(regex, str))
    Species(m[1], Charge(m[2]), m[3:end]...)
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


function Base.show(io::IO, sp::Species)
    out = gas(sp)
    open_bracket = false
    if !(charge(sp) isa Neutral)
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

# Base.show(io::IO, sp::Species) = Base.print(io, sp)
