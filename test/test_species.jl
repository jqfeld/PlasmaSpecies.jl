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
        )
        @test gas(charged) == DiNitrogen()
        @test charge(charged) == Positive(2)
        @test string(electronic_state(charged)) == "B"
        @test vibrational_state(charged) == "1"
        @test isnothing(rotational_state(charged))
    end

    @testset "String parsing" begin
        parsed = Species(" N2(+, X, v=2, J=1) ")
        @test gas(parsed) == DiNitrogen()
        @test charge(parsed) == Positive(1)
        @test string(electronic_state(parsed)) == "X"
        @test vibrational_state(parsed) == "2"
        @test rotational_state(parsed) == "1"

        electron = Species("e")
        @test gas(electron) == Electron()
        @test charge(electron) == Negative(1)

        unknown = Species("Ar(+)")
        @test gas(unknown) isa StringGas
        @test string(gas(unknown)) == "Ar"
        @test charge(unknown) == Positive(1)
    end

    @testset "Parent relationships" begin
        full = Species("N2(+,B,v=3,J=5)")
        vib = Species("N2(+,B,v=3)")
        elec = Species("N2(+,B)")
        base = Species("N2(+)")

        @test get_parent_species(full) == vib
        @test get_parent_species(vib) == elec
        @test get_parent_species(elec) == base
        @test get_parent_species(base) === nothing

        @test is_parent_species(full, vib)
        @test is_parent_species(full, elec)
        @test is_parent_species(full, base)
        @test is_parent_species(vib, full)
        @test is_parent_species(full, Species("N2(B)"))
    end

    @testset "Equality and ordering" begin
        @test Species("N2(X,v=1)") == Species("N2( X , v=1)")

        ordering = Species.(["N2(+)", "N2", "e"])
        @test sort(ordering) == Species.(["N2", "N2(+)", "e"])

        states = Species.(["N2(B)", "N2(X)"])
        @test sort(states) == Species.(["N2(B)", "N2(X)"])
    end

    @testset "Display" begin
        sp = Species("N2(+,X,v=2,J=3)")
        @test string(sp) == "N2(+,X,v=2,J=3)"

        neutral = Species("N2(X,v=2)")
        @test string(neutral) == "N2(X,v=2)"

        bare = Species("N2")
        @test string(bare) == "N2"
    end

    @testset "Mass" begin
        neutral = Species("N2")
        positive = Species("N2(+)")
        negative = Species("N2(-)")
        electron = Species("e")

        me = mass(Electron())
        m_n2 = mass(DiNitrogen())

        @test mass(neutral) == m_n2
        @test mass(positive) ≈ m_n2 - me atol=eps(m_n2)
        @test mass(negative) ≈ m_n2 + me atol=eps(m_n2)
        @test mass(electron) == me
    end
end
