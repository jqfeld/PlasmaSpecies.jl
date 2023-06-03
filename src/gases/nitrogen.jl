"""
    Nitrogen <: Gas

Atomic nitrogen. 
Registered label: `N`

"""
struct Nitrogen <: Gas end
Base.show(io::IO, ::Nitrogen) = print(io,"N")

"""
    DiNitrogen <: Gas

Molecular nitrogen. 
Registered label: `N2`

"""
struct DiNitrogen <: Gas end
Base.show(io::IO, ::DiNitrogen) = print(io,"N2")


"""
    TriNitrogen <: Gas

Tri atomic nitrogen.  
Registered label: `N3`

"""
struct TriNitrogen <: Gas end
Base.show(io::IO, ::TriNitrogen) = print(io,"N3")

"""
    TetraNitrogen <: Gas

Tetra atomic nitrogen.  
Registered label: `N4`

"""
struct TetraNitrogen <: Gas end
Base.show(io::IO, ::TetraNitrogen) = print(io,"N4")

