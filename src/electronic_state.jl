abstract type ElectronicState end

ElectronicState(::Nothing) = nothing
ElectronicState(x::S) where S<: ElectronicState = x

"""
    label(s::ElectronicState) -> Symbol

The conventional letter label of an electronic state (`:X`, `:A`, `:B`, `:a''`, ...).
This is what identifies/matches a state everywhere else in the package (LoKI-B
parsing round-trip, `==`, tree lookup, sorting) — richer physics carried by
[`TermSymbol`](@ref) is additive annotation on top of it, the same way `Species`'s
`energy`/`degeneracy` don't participate in `Species` equality either.
"""
function label end

struct StringElectronicState <: ElectronicState
    name::Symbol
end
StringElectronicState(s::AbstractString) = StringElectronicState(Symbol(s))
label(s::StringElectronicState) = s.name

"""
    TermSymbol <: ElectronicState

An electronic state carrying real term-symbol data alongside its conventional
`label`: spin `multiplicity` (2S+1), orbital angular momentum projection `Λ`
(0=Σ, 1=Π, 2=Δ, 3=Φ, 4=Γ), `parity` (`:g`/`:u`), and `reflection` symmetry
(`:plus`/`:minus`, meaningful for Σ states, and not always available even then).
All but `label` are optional — real data doesn't always carry all of them.
"""
Base.@kwdef struct TermSymbol <: ElectronicState
    label::Symbol
    multiplicity::Union{Int,Nothing} = nothing
    Λ::Union{Int,Nothing} = nothing
    parity::Union{Symbol,Nothing} = nothing
    reflection::Union{Symbol,Nothing} = nothing
end
label(s::TermSymbol) = s.label

const GREEK_LAMBDA_NAMES = ["Sigma" => 0, "Pi" => 1, "Delta" => 2, "Phi" => 3, "Gamma" => 4]

"""
    ElectronicState(x::AbstractString) -> ElectronicState

Parse an electronic-state label. Plain LoKI-B labels (`"X"`, `"A"`, `"a''"`, no
digit suffix) become a [`StringElectronicState`](@ref). Fused term-symbol strings
in the format ExoMol's `ElecState` column uses —
`<label><multiplicity><Λ-name><parity>[<reflection>]`, e.g. `"B3Pig"`,
`"A3Sigmau+"` — decompose into a [`TermSymbol`](@ref).
"""
function ElectronicState(x::AbstractString)
    m = match(r"^([A-Za-z']+?)(\d+)(Sigma|Pi|Delta|Phi|Gamma)([gu])([+-])?$", x)
    isnothing(m) && return StringElectronicState(x)
    label_str, mult, Λname, par, refl = m
    return TermSymbol(
        label = Symbol(label_str),
        multiplicity = parse(Int, mult),
        Λ = Dict(GREEK_LAMBDA_NAMES)[Λname],
        parity = Symbol(par),
        reflection = isnothing(refl) ? nothing : (refl == "+" ? :plus : :minus),
    )
end

Base.:(==)(a::ElectronicState, b::ElectronicState) = label(a) == label(b)
Base.print(io::IO, s::ElectronicState) = print(io, label(s))
Base.show(io::IO, s::ElectronicState) = print(io, label(s))

# order electronic states alphabetically by label
Base.isless(a::ElectronicState, b::ElectronicState) = isless(string(label(a)), string(label(b)))

# concrete electronic state comes before unspecified electronic state
Base.isless(a::ElectronicState, ::Nothing) = false
Base.isless(::Nothing, ::ElectronicState) = true
