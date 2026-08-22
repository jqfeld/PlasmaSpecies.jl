"""
    Electron <: Gas

The electron. The only entry in [`REGISTERED_GASES`](@ref), so the label `"e"`
is resolved before any formula parsing.

A `Species` built on it is forced to `Negative(1)`, and [`mass`](@ref) returns
the electron rest mass directly instead of subtracting it from a neutral.
"""
struct Electron <: Gas end
Base.show(io::IO, ::Electron) = print(io,"e")

"""
    M_E_KG

The electron rest mass in kilograms, CODATA 2018 — the same source as
[`U_KG`](@ref) and the nuclide tables, so an ion mass stays consistent with the
neutral it is derived from.
"""
const M_E_KG = 9.1093837015e-31

mass(::Electron) = M_E_KG
