```@meta; 
CurrentModule = PlasmaSpecies
```
# Species

A `Species` is represented by the following struct
```@docs; canonical=false
Species
```

Usually, it is more convenient to define the species by a string in the format
given defined by [LoKI-B](https://github.com/IST-Lisbon/LoKI). 
For this purpose, `PlasmaSpecies` supplies a convenience constructor: 
```@docs; canonical=false
Species(::String)
```

With this the definition of a species is as easy as
```@example
using PlasmaSpecies
Species("N2(X,v=0,J=10)")
```

A valid string consists of a label for the gas (`N2`) optionally more
detailed information can be given in brackets separated by commas, starting with the charge as a
variable number of `+` or `-` signs or a zero. If no charge is explicitly given
the default (`Neutral()`) is used. The charge is followed by an optional label for the
electronic state. The labels for the vibrational and rotational states are
prefixed by `v=` or `J=`, respectively.
Some valid examples are:
```@repl
using PlasmaSpecies
Species("e")
Species("N3(+)")
Species("N2(+,B)")
Species("N2(-,X,v=2)")
p"N4(+)"
```

The last repl command uses the `p"..."` string macro which can be used to
define species and reactions as well.


## Methods
Various methods exist to query information from a `Species` object.
These include methods to obtain the fields of the struct:
```@repl species-methods; continued = true
using PlasmaSpecies
sp = Species("N2(+,B)")
gas(sp)
charge(sp)
electronic_state(sp)
vibrational_state(sp) === nothing
rotational_state(sp) === nothing
```

Further, the mass of a species (in kg) can be queried if its corresponding `gas` has
the corresponding method defined.`PlasmaSpecies` contains a few registered
gases which are automatically recognized in the species strings and have their
masses defined. The masses of the missing or additional
electrons corresponding to the charge of the species are accounted. 


```@repl species-methods
mass(sp)
mass(p"e")
mass(p"N2")
mass(p"N2") - mass(p"e") ≈ mass(sp)
mass(p"unknown_gas(+,X)")
struct UnknownGas <: Gas end
PlasmaSpecies.mass(::UnknownGas) = π * 1e-26
mass(Species(UnknownGas()))

```
