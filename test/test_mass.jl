using PlasmaSpecies
using Test

const U = 1.66053906660e-27   # dalton in kg

# Checked against reference weights in daltons, not against each other, so a
# wrong constant or table entry cannot cancel out.
@testset "natural-abundance masses" begin
    for (formula, u) in (
        "N"   => 14.007,
        "O"   => 15.999,
        "H"   => 1.008,
        "C"   => 12.011,
        "N2"  => 28.014,
        "O2"  => 31.998,
        "O3"  => 47.997,
        "Ar"  => 39.95,
        "He"  => 4.002602,
        "H2O" => 18.015,
        "CO2" => 44.009,
        "CO"  => 28.010,
        "CH4" => 16.043,
        "OH"  => 17.007,
        "NO"  => 30.006,
        "NH3" => 17.031,
        "SF6" => 146.05,
    )
        @test mass(Species(formula)) ≈ u * U rtol = 1.0e-4
    end
end

@testset "electron mass" begin
    # CODATA 2018, the same source as U and the nuclide tables.
    @test mass(Species("e")) ≈ 9.1093837015e-31 rtol = 1.0e-12

    # Independent cross-check: the electron's relative atomic mass in daltons,
    # converted with U, rather than the kg value compared against itself.
    @test mass(Species("e")) ≈ 5.48579909065e-4 * U rtol = 1.0e-10
end

@testset "structure and charge" begin
    @test mass(Species("N2")) ≈ 2 * mass(Species("N"))
    @test mass(Species("N2[X,vib=0]")) == mass(Species("N2"))
    @test mass(Species("N2[X,vib=0]")) == mass(p"N2")

    # An ion differs from its neutral by exactly one electron mass per charge.
    @test mass(Species("N2[-]")) ≈ mass(Species("N2")) + mass(Species("e"))
    @test mass(Species("N2[+]")) ≈ mass(Species("N2")) - mass(Species("e"))
    @test mass(Species("O2[-]")) ≈ mass(Species("O2")) + mass(Species("e"))
    @test mass(Species("O2[++]")) ≈ mass(Species("O2")) - 2mass(Species("e"))
end

@testset "isotopes" begin
    @test mass(Species("14N")) ≈ 14.00307400443 * U rtol = 1.0e-9
    @test mass(Species("12C")) ≈ 12.0 * U rtol = 1.0e-12
    @test mass(Species("14N2")) ≈ 2 * 14.00307400443 * U rtol = 1.0e-9

    # Lighter than the natural mixture, which carries heavier minor isotopes.
    @test mass(Species("14N2")) < mass(Species("N2"))

    # Deuterium and tritium resolve without a mass number.
    @test mass(Species("D")) ≈ 2.01410177784 * U rtol = 1.0e-9
    @test mass(Species("T")) ≈ 3.01604927790 * U rtol = 1.0e-9
    @test mass(Species("D2")) ≈ 2 * mass(Species("D"))

    # Independent cross-check: ExoMol's own .def file gives SiO's isotopologue
    # mass as 43.971842 Da (asserted in SpectraUtils.jl's test_exomol.jl).
    @test mass(Species("28Si-16O")) ≈ 43.971842 * U rtol = 1.0e-6
    @test mass(Species("16O-1H")) ≈ (15.99491461957 + 1.00782503190) * U rtol = 1.0e-9
    @test mass(Species("12C-16O2")) ≈ (12.0 + 2 * 15.99491461957) * U rtol = 1.0e-9
end

@testset "unresolvable gases still throw" begin
    @test gas(Species("Xx")) isa StringGas
    @test_throws ErrorException mass(Species("Xx"))
    # Unknown *isotope* of a known element: the element table has N, the
    # isotope table has no N-99.
    @test gas(Species("99N")) isa StringGas
    @test_throws ErrorException mass(Species("99N"))
end
