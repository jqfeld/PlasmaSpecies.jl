"""
    Nuclide

A composition primitive: either an [`Element`](@ref) at natural isotopic
abundance or a specific [`Isotope`](@ref). Gases are compositions over
nuclides — see [`Molecule`](@ref).
"""
abstract type Nuclide end

"""
    U_KG

The unified atomic mass unit (dalton) in kilograms, CODATA 2018.
"""
const U_KG = 1.66053906660e-27

"""
    Element(symbol, Z, weight) <: Nuclide

A chemical element at natural isotopic abundance. `weight` is the standard
atomic weight in daltons (IUPAC 2021 conventional values).

A bare element symbol in a formula resolves to this, so `"N2"` is two
natural-abundance nitrogens (2 × 14.007 u), which is what LoKI-B and BOLSIG+
assume for a bulk gas. Give a mass number (`"14N2"`) for a single nuclide.
"""
struct Element <: Nuclide
  symbol::Symbol
  Z::Int
  weight::Float64
end

"""
    Isotope(symbol, Z, A, mass) <: Nuclide

A single nuclide. `mass` is the exact atomic mass in daltons (AME2020).
"""
struct Isotope <: Nuclide
  symbol::Symbol
  Z::Int
  A::Int
  mass::Float64
end

"""
    mass(n::Nuclide) -> Float64

Mass of one atom in kilograms.
"""
mass(n::Element) = n.weight * U_KG
mass(n::Isotope) = n.mass * U_KG

"""
    symbol(n::Nuclide) -> Symbol

Element symbol, with no mass number attached (`:O` for both `O` and `¹⁶O`).
"""
symbol(n::Element) = n.symbol
symbol(n::Isotope) = n.symbol

"""
    atomic_number(n::Nuclide) -> Int
"""
atomic_number(n::Nuclide) = n.Z

"""
    mass_number(n::Nuclide) -> Union{Int,Nothing}

Nucleon count for an [`Isotope`](@ref), `nothing` for an [`Element`](@ref).
"""
mass_number(::Element) = nothing
mass_number(n::Isotope) = n.A

Base.print(io::IO, n::Element) = print(io, n.symbol)
Base.print(io::IO, n::Isotope) = print(io, n.A, n.symbol)
Base.show(io::IO, n::Element) = print(io, n.symbol)
Base.show(io::IO, n::Isotope) = print(io, superscript(n.A), n.symbol)

# By (Z, mass number), so canonical compositions sort reproducibly.
Base.isless(a::Nuclide, b::Nuclide) =
  isless((a.Z, something(mass_number(a), 0)), (b.Z, something(mass_number(b), 0)))

const SUPERSCRIPTS = Dict(collect("0123456789") .=> collect("⁰¹²³⁴⁵⁶⁷⁸⁹"))
const SUBSCRIPTS = Dict(collect("0123456789") .=> collect("₀₁₂₃₄₅₆₇₈₉"))

superscript(n::Integer) = map(c -> SUPERSCRIPTS[c], string(n))
subscript(n::Integer) = map(c -> SUBSCRIPTS[c], string(n))
