abstract type Charge end
# NOTE: I am not sure the separation of Positive, Negative and Neutral is
# necessarry. Probably, Charge{Val} would be enough, but maybe with the three
# types some nice dispatch magic might be possible later.
struct Positive{Val} <: Charge end
Positive(c) = Positive{Val{c}}()
struct Negative{Val} <: Charge end
Negative(c) = Negative{Val{c}}()
struct Neutral <: Charge end
abstract type Species{C<:Charge} end
abstract type ElectronicState end
abstract type VibrationalState end
abstract type RotationalState end
struct RoVibrationalState{V<:VibrationalState,R<:RotationalState}
    vibrational_state::V
    rotational_state::R
end
rovibrational_state(::S) where S<:Species = RoVibrationalState{vibrational_state(S()), rotational_state(S())}

charge_type(::S) where {C,S<:Species{C}} = C

charge(::Positive{Val{C}}) where C = C
name(::Positive{Val{C}}) where C = '+'^C
charge(::Negative{Val{C}}) where C = -C
name(::Negative{Val{C}}) where C = '-'^C
charge(::Neutral) = 0
name(::Neutral) = "0"
charge(::Species{C}) where C = charge(C)

name(::S) where {C,S<:Species{C}} = name(S) * "($(name(C())))"
Base.show(io::IO, s::Species) = show(io,name(s))

ground_state(s::Species) = error("No ground state defined for $(name(s)).")

energy(::S, ::E, rv::RoVibrationalState{V,R}=zero(rotational_state(S()))) where 
    {S<:Species, E<:ElectronicState, V<:VibrationalState, R<:RotationalState} = 
    error("energy() not implemented for this species or state")

energy(::S,rv::RoVibrationalState{V,R}=zero(rotational_state(S()))) where 
    {S<:Species,  V<:VibrationalState, R<:RotationalState} = 
    energy(S(), ground_state(S()), rv)



spin(::ElectronicState) = error("spin() needs to be implemented for an electronic state")
total_angular_momentum(::ElectronicState) = error("total_angular_momentum() needs to be implemented for an electronic state")
orbital_angular_momentum(::ElectronicState) = error("orbital_angular_momentum() needs to be implemented for an electronic state")


include("species/Nitrogen.jl")
