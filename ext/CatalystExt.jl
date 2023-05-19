module CatalystExt

using PlasmaSpecies
using Catalyst
function PlasmaSpecies.to_catalyst(sp::Species)
    symbol = Symbol(string(sp))
    return (@species $symbol)[1]
end

end
