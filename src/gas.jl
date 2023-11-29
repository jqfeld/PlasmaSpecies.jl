abstract type Gas end


include("gases/electron.jl")
include("gases/nitrogen.jl")
const REGISTERED_GASES = Dict(
    "e" => Electron,
    "N" => Nitrogen,
    "N2" => DiNitrogen,
    "N3" => TriNitrogen,
    "N4" => TetraNitrogen,
)
Gas(x) = haskey(REGISTERED_GASES, x) ? REGISTERED_GASES[x]() : StringGas(x)

struct StringGas <: Gas
    name::String
end

Base.show(io::IO, g::StringGas) = show(io, g.name)
Base.print(io::IO, g::StringGas) = print(io, g.name)

Base.isless(a::Gas, b::Gas) = isless(string(a), string(b))

mass(::Gas) = error(
    """Mass needs to be defined. If you are using a p\"\" string, the gas seems to be unknown to plasma species. 
    In this case either post an issue in the git repo or define a new Gas struct with corresponding mass(::Gas) function.""")
