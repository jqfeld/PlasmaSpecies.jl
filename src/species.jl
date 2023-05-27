
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

Base.print(io::IO, c::Positive) = Base.print(io::IO, '+'^c.value)
Base.print(io::IO, c::Negative) = Base.print(io::IO, '-'^c.value)
Base.print(io::IO, ::Neutral) = Base.print(io::IO, "")


struct Species{S,C<:Charge,E,V,J}
    gas::S
    charge::C
    electronic_state::E
    vibrational_state::V
    rotational_state::J
end
isnegative(x) = false
isnegative(::Species{S,Negative,E,V,J}) where {S,E,V,J} = true
ispositive(x) = false
ispositive(::Species{S,Positive,E,V,J}) where {S,E,V,J} = true
isneutral(x) = false
isneutral(::Species{S,Neutral,E,V,J}) where {S,E,V,J} = true

Base.:(==)(x::S, y::S) where {S<:Species} =
    gas(x) == gas(y) &&
    charge(x) == charge(y) &&
    electronic_state(x) == electronic_state(y) &&
    vibrational_state(x) == vibrational_state(y) &&
    rotational_state(x) == rotational_state(y)

function Species(str::String)
    regex = r"(^(?:\w|\d)*)(?:\(([+\-]+)?,?([^,\)]+)?(?:,v=([^,\)]*)(?:,J=([^,\)]*))?)?\))?"
    m = collect(match(regex, str))
    Species(m[1], m[1] == "e" ? Negative(1) : Charge(m[2]), m[3:end]...)
end


gas(sp::Species) = sp.gas
charge(sp::Species) = sp.charge
electronic_state(sp::Species) = sp.electronic_state
vibrational_state(sp::Species) = sp.vibrational_state
rotational_state(sp::Species) = sp.rotational_state

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
    else
        return nothing
    end
end

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


function Base.show(io::IO, sp::Species)
    out = gas(sp)
    open_bracket = false
    if !(charge(sp) isa Neutral) && gas(sp) != "e"
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

# Sorting
Base.isless(::Negative, ::Charge )= false

Base.isless(::Positive, ::Negative) = true
Base.isless(::Neutral, ::Negative) = true

Base.isless(::Positive, ::Neutral) = false
Base.isless(::Neutral, ::Positive) = true

Base.isless(a::Negative, b::Negative) = isless(a.value, b.value)
Base.isless(a::Positive, b::Positive) = isless(a.value, b.value)

function Base.isless(a::Species, b::Species)
    if gas(a) == "e" 
        return false
    elseif gas(b) == "e"
        return true
    end

    if charge(a) != charge(b)
        return isless(charge(a), charge(b))
    end
    
    if gas(a) != gas(b)
        return isless(gas(a), gas(b))
    end

    if electronic_state(a) != electronic_state(b)
        if electronic_state(a) == "X"
            return true
        end
        if electronic_state(b) == "X"
            return false
        end
        return isless(electronic_state(a), electronic_state(b))
    end

    if vibrational_state(a) != vibrational_state(b)
        if isnothing(vibrational_state(a))
            return false
        end
        if isnothing(vibrational_state(b))
            return true
        end
        return isless(parse(Int,vibrational_state(a)),parse(Int, vibrational_state(b)))
    end

    if rotational_state(a) != rotational_state(b)
        if isnothing(rotational_state(a))
            return false
        end
        if isnothing(rotational_state(b))
            return true
        end
        return isless(parse(Int,rotational_state(a)), parse(Int,rotational_state(b)))
    end
end


