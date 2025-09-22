"""
    Oxygen <: Gas

Atomic oxygen.
Registered label: `O`

"""
struct Oxygen <: Gas end
Base.show(io::IO, ::Oxygen) = print(io,"O")
Base.print(io::IO, ::Oxygen) = print(io, "O")
mass(::Oxygen) = 16*1.67e-27 # kg

"""
    DiOxygen <: Gas

Molecular oxygen.
Registered label: `O2`

"""
struct DiOxygen <: Gas end
Base.show(io::IO, ::DiOxygen) = print(io,"O₂")
Base.print(io::IO, ::DiOxygen) = print(io,"O2")
mass(::DiOxygen) = mass(Oxygen())*2 # kg


"""
    TriOxygen <: Gas

Tri atomic oxygen (ozone).
Registered label: `O3`

"""
struct TriOxygen <: Gas end
Base.show(io::IO, ::TriOxygen) = print(io,"O₃")
Base.print(io::IO, ::TriOxygen) = print(io, "O3")
mass(::TriOxygen) = mass(Oxygen())*3 # kg

"""
    TetraOxygen <: Gas

Tetra atomic oxygen.
Registered label: `O4`

"""
struct TetraOxygen <: Gas end
Base.show(io::IO, ::TetraOxygen) = print(io,"O₄")
Base.print(io::IO, ::TetraOxygen) = print(io, "O4")
mass(::TetraOxygen) = mass(Oxygen())*4 # kg

