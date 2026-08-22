using PlasmaSpecies
using Test

@testset "ElectronicState" begin
    @testset "Plain labels stay StringElectronicState" begin
        for s in ["X", "A", "B", "a", "a'", "a''", "w", "W"]
            es = ElectronicState(s)
            @test es isa StringElectronicState
            @test label(es) == Symbol(s)
            @test string(es) == s
        end
    end

    @testset "Fused term-symbol strings (real ExoMol N2 ElecState values)" begin
        # From ~/.julia/scratchspaces/.../14N2__WCCRMT.states.bz2 (this session's
        # cached ExoMol download) — every distinct ElecState value in the real file.
        cases = [
            ("A3Sigmau+", :A, 3, 0, :u, :plus),
            ("B'3Sigmau", Symbol("B'"), 3, 0, :u, nothing),
            ("B3Pig",     :B,  3, 1, :g, nothing),
            ("C3Piu",     :C,  3, 1, :u, nothing),
            ("W3Deltau",  :W,  3, 2, :u, nothing),
        ]
        for (str, lbl, mult, Λ, par, refl) in cases
            es = ElectronicState(str)
            @test es isa TermSymbol
            @test label(es) == lbl
            @test es.multiplicity == mult
            @test es.Λ == Λ
            @test es.parity == par
            @test es.reflection == refl
        end
    end

    @testset "Round-trip identity is label-only" begin
        # print/show/string stay label-only for LoKI-B round-trip — the term-symbol
        # fields are annotation, not part of the species' textual identity.
        es = ElectronicState("A3Sigmau+")
        @test string(es) == "A"
        @test sprint(show, es) == "A"
    end

    @testset "Cross-type equality keys on label" begin
        @test TermSymbol(label=:B, multiplicity=3, Λ=1, parity=:g) == StringElectronicState("B")
        @test StringElectronicState("B") == TermSymbol(label=:B)
        @test TermSymbol(label=:B) != StringElectronicState("C")
    end

    @testset "Species round-trip unaffected by the parser upgrade" begin
        sp = Species("N2[X,vib=2,rot=3]")
        @test electronic_state(sp) isa StringElectronicState
        @test string(sp) == "N2[X,vib=2,rot=3]"
    end

    @testset "Sorting still alphabetical by label" begin
        states = ElectronicState.(["B", "A3Sigmau+", "C3Piu"])
        @test label.(sort(states)) == [:A, :B, :C]
    end
end
