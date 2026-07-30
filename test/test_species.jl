using PlasmaSpecies
using PlasmaSpecies: Electron
using Test

@testset "Species" begin
    @testset "Constructors" begin
        neutral = Species(DiNitrogen())
        @test gas(neutral) == DiNitrogen()
        @test charge(neutral) isa Neutral
        @test isnothing(electronic_state(neutral))
        @test isnothing(vibrational_state(neutral))
        @test isnothing(rotational_state(neutral))

        charged = Species(
            DiNitrogen(),
            Positive(2),
            StringElectronicState("B"),
            "1",
            nothing,
            nothing,
            nothing,
        )
        @test gas(charged) == DiNitrogen()
        @test charge(charged) == Positive(2)
        @test string(electronic_state(charged)) == "B"
        @test vibrational_state(charged) == "1"
        @test isnothing(rotational_state(charged))
    end

    @testset "Keyword constructor" begin
        sp = Species(gas=DiNitrogen())
        @test gas(sp) == DiNitrogen()
        @test charge(sp) isa Neutral
        @test isnothing(electronic_state(sp))
        @test isnothing(vibrational_state(sp))
        @test isnothing(rotational_state(sp))

        sp2 = Species(gas=DiNitrogen(), charge=Positive(1), electronic_state=StringElectronicState("B"), vibrational_state="3")
        @test gas(sp2) == DiNitrogen()
        @test charge(sp2) == Positive(1)
        @test string(electronic_state(sp2)) == "B"
        @test vibrational_state(sp2) == "3"
        @test isnothing(rotational_state(sp2))
    end

    @testset "String parsing" begin
        parsed = Species(" N2[+, X, vib=2, rot=1] ")
        @test gas(parsed) == DiNitrogen()
        @test charge(parsed) == Positive(1)
        @test string(electronic_state(parsed)) == "X"
        @test vibrational_state(parsed) === 2
        @test rotational_state(parsed) === 1

        electron = Species("e")
        @test gas(electron) == Electron()
        @test charge(electron) == Negative(1)

        unknown = Species("Ar[+]")
        @test gas(unknown) isa StringGas
        @test string(gas(unknown)) == "Ar"
        @test charge(unknown) == Positive(1)
    end

    @testset "Quantum label parsing" begin
        @test vibrational_state(Species("N2[X,vib=3]")) === 3
        @test vibrational_state(Species("N2[X,vib=0-4]")) === 0:4
        @test vibrational_state(Species("N2[X,vib=(1,2)]")) === (1, 2)
        @test vibrational_state(Species("N2[X,vib=(v1=1,v2=0)]")) === (v1=1, v2=0)
        @test vibrational_state(Species("N2[X,vib=combined]")) === :combined

        @test rotational_state(Species("N2[X,vib=0,rot=5]")) === 5
        @test rotational_state(Species("N2[X,vib=0,rot=(1,2)]")) === (1, 2)
    end

    @testset "Parent relationships" begin
        full = Species("N2[+,B,vib=3,rot=5]")
        vib = Species("N2[+,B,vib=3]")
        elec = Species("N2[+,B]")
        base = Species("N2[+]")

        @test get_parent_species(full) == vib
        @test get_parent_species(vib) == elec
        @test get_parent_species(elec) == base
        @test get_parent_species(base) === nothing

        @test is_parent_species(full, vib)
        @test is_parent_species(full, elec)
        @test is_parent_species(full, base)
        @test is_parent_species(vib, full)
        @test is_parent_species(full, Species("N2[B]"))

        ranged = Species("N2[+,B,vib=0-4]")
        @test is_parent_species(vib, ranged)
        @test !is_parent_species(Species("N2[+,B,vib=5]"), ranged)
        @test !is_parent_species(Species("N2[+,A,vib=3]"), ranged)
    end

    @testset "Base.in" begin
        full   = Species("N2[+,B,vib=3,rot=5]")
        vib    = Species("N2[+,B,vib=3]")
        elec   = Species("N2[+,B]")
        base   = Species("N2[+]")
        ranged = Species("N2[+,B,vib=0-4]")

        # reflexive, unlike is_parent_species
        @test full in full
        @test !is_parent_species(full, full)

        # nothing fields on the right-hand side are wildcards
        @test full in vib
        @test full in elec
        @test full in base
        @test vib in elec

        # a mismatched concrete field is not a wildcard
        @test !(full in Species("N2[+,B,vib=2,rot=5]"))
        @test !(full in Species("N2[+,A]"))

        # range containment, and rotational_state actually checked (unlike is_parent_species)
        @test vib in ranged
        @test !(Species("N2[+,B,vib=5]") in ranged)
        @test full in Species("N2[+,B,vib=3,rot=5]")
        @test !(Species("N2[+,B,vib=3,rot=6]") in Species("N2[+,B,vib=3,rot=5]"))
    end

    @testset "Equality and ordering" begin
        @test Species("N2[X,vib=1]") == Species("N2[ X , vib=1]")

        ordering = Species.(["N2[+]", "N2", "e"])
        @test sort(ordering) == Species.(["N2", "N2[+]", "e"])

        states = Species.(["N2[B]", "N2[X]"])
        @test sort(states) == Species.(["N2[B]", "N2[X]"])

        # vibrational_state is now a real Int (not a string needing `parse`), so
        # comparing species that differ only by vib level must not throw.
        vibs = Species.(["N2[X,vib=2]", "N2[X,vib=0]", "N2[X,vib=1]"])
        @test sort(vibs) == Species.(["N2[X,vib=0]", "N2[X,vib=1]", "N2[X,vib=2]"])
    end

    @testset "Display" begin
        sp = Species("N2[+,X,vib=2,rot=3]")
        @test string(sp) == "N2[+,X,vib=2,rot=3]"

        neutral = Species("N2[X,vib=2]")
        @test string(neutral) == "N2[X,vib=2]"

        bare = Species("N2")
        @test string(bare) == "N2"

        ranged = Species("N2[X,vib=0-4]")
        @test string(ranged) == "N2[X,vib=0-4]"
        @test Species(string(ranged)) == ranged
    end

    @testset "Mass" begin
        neutral = Species("N2")
        positive = Species("N2[+]")
        negative = Species("N2[-]")
        electron = Species("e")

        me = mass(Electron())
        m_n2 = mass(DiNitrogen())

        @test mass(neutral) == m_n2
        @test mass(positive) ≈ m_n2 - me atol=eps(m_n2)
        @test mass(negative) ≈ m_n2 + me atol=eps(m_n2)
        @test mass(electron) == me
    end
end
