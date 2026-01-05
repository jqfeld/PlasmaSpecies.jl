abstract type ElectronicState end

ElectronicState(::Nothing) = nothing
ElectronicState(x::S) where S<: ElectronicState = x
ElectronicState(x::AbstractString) = StringElectronicState(x)

struct StringElectronicState <: ElectronicState
    name::String
end
Base.show(io::IO, g::StringElectronicState) = show(io, g.name)
Base.print(io::IO, g::StringElectronicState) = print(io, g.name)

# order electronic states alphabetically
Base.isless(a::ElectronicState, b::ElectronicState) = isless(string(a), string(b))

# concrete electronic state comes before unspecified electronic state
Base.isless(a::ElectronicState, ::Nothing) = false
Base.isless(::Nothing, ::ElectronicState) = true
