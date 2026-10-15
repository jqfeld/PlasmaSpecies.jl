using PlasmaSpecies
using Test

@testset "plain formulas" begin
    for (s, expected) in (
        "Ar"  => [(:Ar, 1)],
        "N2"  => [(:N, 2)],
        "CO2" => [(:C, 1), (:O, 2)],
        "H2O" => [(:H, 2), (:O, 1)],
        "CH4" => [(:C, 1), (:H, 4)],
        "OH"  => [(:O, 1), (:H, 1)],
        "NO"  => [(:N, 1), (:O, 1)],
        "SF6" => [(:S, 1), (:F, 6)],
        "N2O" => [(:N, 2), (:O, 1)],
    )
        g = parse_formula(s)
        @test g isa Molecule
        # written order is preserved
        @test [(PlasmaSpecies.symbol(n), c) for (n, c) in composition(g)] == expected
        @test string(g) == s
    end
end

@testset "greedy two-letter symbols" begin
    # Case decides: CO is carbon monoxide, Co is cobalt.
    symbols(s) = collect(PlasmaSpecies.symbol.(first.(composition(parse_formula(s)))))
    @test symbols("CO") == [:C, :O]
    @test symbols("Co") == [:Co]
    @test symbols("NO") == [:N, :O]
    @test symbols("No") == [:No]
    # Ne is neon, not N + e.
    @test symbols("Ne") == [:Ne]
end

@testset "braced mass numbers" begin
    for s in ("{14}N2", "{16}O{1}H", "{12}C{16}O2", "{28}Si{16}O", "{13}C{1}H4", "{16}O", "{13}CO2")
        g = parse_formula(s)
        @test g isa Molecule
        @test string(g) == s          # round-trips in the braced spelling
        @test isfinite(mass(g))
    end
    # A braced number is the mass number of the symbol after it, a trailing
    # digit is a count.
    g = parse_formula("{12}C{16}O2")
    @test [(PlasmaSpecies.symbol(n), PlasmaSpecies.mass_number(n), c) for (n, c) in composition(g)] ==
          [(:C, 12, 1), (:O, 16, 2)]
    # Natural-abundance nuclides have no mass number, and mix with labelled ones.
    @test all(isnothing ∘ PlasmaSpecies.mass_number ∘ first, composition(parse_formula("CO2")))
    @test [PlasmaSpecies.mass_number(n) for (n, _) in composition(parse_formula("{13}CO2"))] == [13, nothing]
    # A brace must label a symbol.
    for s in ("{14}", "{14}2", "{}N", "N{14}", "{a}N")
        @test isnothing(parse_formula(s))
    end
end

@testset "legacy isotope spelling is rejected" begin
    for s in ("14N2", "16O-1H", "12C-16O2", "O-1H", "16O")
        @test isnothing(parse_formula(s))
        @test_throws ArgumentError Gas(s)
        @test_throws ArgumentError Species(s)
    end
    # Labels that merely contain digits or dashes are still StringGas.
    for s in ("N2-", "Ar-site", "x1", "-")
        @test Gas(s) isa StringGas
    end
end

@testset "ExoMol iso_slug interop" begin
    for (slug, braced) in ("16O-1H" => "{16}O{1}H", "12C-16O2" => "{12}C{16}O2",
                           "28Si-16O" => "{28}Si{16}O", "14N2" => "{14}N2")
        g = from_iso_slug(slug)
        @test g == parse_formula(braced)
        @test string(g) == braced
        @test iso_slug(g) == slug
    end
    @test_throws ArgumentError from_iso_slug("O-1H")      # no mass number
    @test_throws ArgumentError from_iso_slug("99N")       # unknown nuclide
    @test_throws ArgumentError iso_slug(parse_formula("{13}CO2"))  # not fully resolved
end

@testset "fallback to StringGas" begin
    # unknown element, unknown isotope, malformed, non-formula decoration
    for s in ("Xx", "99N", "Ar*", "N2-", "2", "-", "n2", "")
        @test isnothing(parse_formula(s))
    end
    @test Gas("Xx") isa StringGas
    @test Gas("N2-") isa StringGas
    # the electron keeps its special, non-formula label
    @test Gas("e") isa Electron
end

@testset "equality is canonical, display is as written" begin
    # Same composition written two ways is one gas, and one Dict key.
    @test Species("H2O") == Species("OH2")
    @test hash(Species("H2O")) == hash(Species("OH2"))
    @test haskey(Dict(Species("H2O") => 1), Species("OH2"))
    # but each still prints the way it was written
    @test string(gas(Species("H2O"))) == "H2O"
    @test string(gas(Species("OH2"))) == "OH2"

    # Repeated groups merge in the canonical form: CH3CH3 is C2H6.
    @test Species("CH3CH3") == Species("C2H6")
    @test mass(Species("CH3CH3")) == mass(Species("C2H6"))

    # Different gases stay different.
    @test Species("CO") != Species("Co")
    @test Species("N2") != Species("{14}N2")
end

@testset "pretty display" begin
    @test sprint(show, gas(Species("N2"))) == "N₂"
    @test sprint(show, gas(Species("CO2"))) == "CO₂"
    @test sprint(show, gas(Species("{16}O{1}H"))) == "¹⁶O¹H"
    # print stays plain ASCII — LoKI-B round-tripping depends on it
    @test sprint(print, gas(Species("N2"))) == "N2"
    @test sprint(print, gas(Species("{16}O{1}H"))) == "{16}O{1}H"
end

@testset "braced mass numbers survive the species parser" begin
    # Braces in the gas token must not cut off the charge or state fields.
    sp = Species("{16}O{1}H[+]")
    @test string(gas(sp)) == "{16}O{1}H"
    @test charge(sp) == Positive(1)

    sp2 = Species("{12}C{16}O2[X,vib=(0,1,0)]")
    @test string(gas(sp2)) == "{12}C{16}O2"
    @test label(electronic_state(sp2)) == :X
    @test vibrational_state(sp2) == (0, 1, 0)

    for s in ("N2", "N2[+]", "N2[+,B,vib=3,rot=1]", "O2[X,vib=0]",
              "CO2[X,vib=(1,0,0)]", "Ar[++]", "H2O", "{16}O{1}H", "{14}N2[+]", "e")
        @test Species(string(Species(s))) == Species(s)
    end
end

@testset "common plasma gases resolve without registration" begin
    for s in ("Ar", "He", "Ne", "Kr", "Xe", "H2", "H2O", "OH", "NO", "N2O",
              "CO", "CO2", "CH4", "NH3", "SF6", "CF4", "Cl2", "SiH4", "C2H2",
              "C2H6", "SO2", "HCl", "F2", "O3", "N4")
        g = Gas(s)
        @test g isa Molecule
        @test isfinite(mass(g))
        @test string(g) == s
    end
end
