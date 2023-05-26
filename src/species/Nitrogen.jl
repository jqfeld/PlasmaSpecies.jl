struct Nitrogen{C} <: Species{C} end
name(::Type{Nitrogen{C}}) where C = "N"


struct DiNitrogen{C} <: Species{C} end
name(::Type{DiNitrogen{C}}) where C = "N2"

DiNitrogen(charge::C=Neutral()) where C = DiNitrogen{C}()

struct DiNitrogenElecState{Val} <: ElectronicState end
ground_state(::DiNitrogen{Neutral}) = DiNitrogenElecState{Val{:X}}

struct DiNitrogenVibState <: VibrationalState
    v::UInt8
end

struct DiNitrogenRotState <: RotationalState
    R::UInt8
end

vibrational_state(::DiNitrogen) = DiNitrogenVibState
rotational_state(::DiNitrogen) = DiNitrogenRotState
