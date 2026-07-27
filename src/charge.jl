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

# `print`/`string` fall back to `show` automatically (no separate `print` needed
# unless a type wants a different plain-vs-pretty form — see DiNitrogen etc.)
Base.show(io::IO, c::Positive) = print(io, '+'^c.value)
Base.show(io::IO, c::Negative) = print(io, '-'^c.value)
Base.show(io::IO, ::Neutral) = print(io, "")

to_value(c::Positive) = c.value
to_value(c::Negative) = -c.value
to_value(c::Neutral) = 0
