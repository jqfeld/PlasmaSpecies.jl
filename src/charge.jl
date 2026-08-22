# NOTE: I am not sure the separation of Positive, Negative and Neutral is
# necessarry. Probably, Charge{Val} would be enough, but maybe with the three
# types some nice dispatch magic might be possible later.
"""
    Charge

Charge state of a [`Species`](@ref): [`Positive`](@ref), [`Negative`](@ref) or
[`Neutral`](@ref). Separate types rather than a signed integer, so charge can be
dispatched on.

    Charge(s::AbstractString) -> Charge
    Charge(::Nothing) -> Neutral()

Parse the charge field of a LoKI-B species string: a run of `+` or `-` gives
`Positive(n)`/`Negative(n)`, an empty string or `"0"` gives `Neutral()`.
Anything else throws.
"""
abstract type Charge end

"""
    Positive(n) <: Charge

`n`-fold positively charged, written as `n` `+` signs (`"N2[+]"`, `"N2[++]"`).
"""
struct Positive <: Charge
    value::Int8
end

"""
    Negative(n) <: Charge

`n`-fold negatively charged, written as `n` `-` signs (`"O[-]"`). An
[`Electron`](@ref) species is always `Negative(1)`.
"""
struct Negative <: Charge
    value::Int8
end

"""
    Neutral() <: Charge

Uncharged. The default when a species string carries no charge field, and
printed as the empty string.
"""
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
# unless a type wants a different plain-vs-pretty form — see Molecule etc.)
Base.show(io::IO, c::Positive) = print(io, '+'^c.value)
Base.show(io::IO, c::Negative) = print(io, '-'^c.value)
Base.show(io::IO, ::Neutral) = print(io, "")

"""
    to_value(c::Charge) -> Int

Signed charge number: `+n`, `-n` or `0`. Used by
[`ischargebalanced`](@ref) and by [`mass`](@ref) to count missing/extra
electrons. Not exported.
"""
to_value(c::Positive) = c.value
to_value(c::Negative) = -c.value
to_value(c::Neutral) = 0
