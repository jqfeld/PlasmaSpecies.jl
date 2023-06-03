abstract type Gas end


include("gases/nitrogen.jl")
const REGISTERED_GASES = Dict(
    "N"   => Nitrogen,
    "N2"  => DiNitrogen,
    "N3"  => TriNitrogen,
    "N4"  => TetraNitrogen,
)
Gas(x) = haskey(REGISTERED_GASES, x) ? REGISTERED_GASES[x]() : StringGas(x)

struct StringGas
    name::String
end

Base.show(io::IO, g::StringGas) = show(io,g.name)
Base.print(io::IO, g::StringGas) = print(io,g.name)

Base.isless(a::Gas, b::Gas) = isless(string(a), string(b))
