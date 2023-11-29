"""
    Electron <: Gas

Atomic nitrogen. 
Registered label: `e`

"""
struct Electron <: Gas end
Base.show(io::IO, ::Electron) = print(io,"e")
mass(::Electron) = 9.11e-31 # kg
