```@meta; 
CurrentModule = PlasmaSpecies
```

# PlasmaSpecies

[PlasmaSpecies.jl](https://github.com/jqfeld/PlasmaSpecies.jl) implements some
basic types for handling the various species in a plasma species. 
Some of its features are:

* Data structures to define species in the plasma (e.g., electrons, ions, and neutrals)
* Methods to construct a tree data structure which shows the relations between
  the various species.
* Methods to define reactions between the species and a way to produce the
  corresponding [Catalyst.jl](https://github.com/SciML/MethodOfLines.jl)
  `Reaction`s.
* Some convenience functionality, e.g., a string macro to define species
  (`p"e"` or `p"N2"`) and reactions (`p"e + N2 --> 2e + N2(+)`) using the
  notation introduced by the [LoKI-B](https://github.com/IST-Lisbon/LoKI)
  Boltzmann solver.

