"""
    Nitrogen <: Gas

Atomic nitrogen. 
Registered label: `N`

"""
struct Nitrogen <: Gas end
Base.show(io::IO, ::Nitrogen) = print(io,"N")
mass(::Nitrogen) = 14*1.67e-27 # kg

"""
    DiNitrogen <: Gas

Molecular nitrogen. 
Registered label: `N2`

"""
struct DiNitrogen <: Gas end
Base.show(io::IO, ::DiNitrogen) = print(io,"N₂")
mass(::DiNitrogen) = mass(Nitrogen())*2 # kg


"""
    TriNitrogen <: Gas

Tri atomic nitrogen.  
Registered label: `N3`

"""
struct TriNitrogen <: Gas end
Base.show(io::IO, ::TriNitrogen) = print(io,"N₃")
mass(::TriNitrogen) = mass(Nitrogen())*3 # kg

"""
    TetraNitrogen <: Gas

Tetra atomic nitrogen.  
Registered label: `N4`

"""
struct TetraNitrogen <: Gas end
Base.show(io::IO, ::TetraNitrogen) = print(io,"N₄")
mass(::TetraNitrogen) = mass(Nitrogen())*4 # kg

