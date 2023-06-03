abstract type ElectronicState end

ElectronicState(::Nothing) = nothing
ElectronicState(x::S) where S<: ElectronicState = x
ElectronicState(x::AbstractString) = StringElectronicState(x)

struct StringElectronicState
    name::String
end
Base.show(io::IO, g::StringElectronicState) = show(io, g.name)
Base.print(io::IO, g::StringElectronicState) = print(io, g.name)

