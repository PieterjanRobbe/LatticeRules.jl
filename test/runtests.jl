using LatticeRules, SpecialFunctions, Statistics, Test

@testset "LatticeRules" begin

    @testset "Constructing LatticeRule" begin

        # test constructor with generating vector, number of dimensions and max number of points
        @testset "LatticeRule(z, s, n)" begin
            lattice_rule = LatticeRule([UInt32(1), UInt32(5)], 2, 8)
            @test ndims(lattice_rule) == 2
            @test length(lattice_rule) == 8
            @test size(lattice_rule) == (8,)
        end

        # test constructor with generating vector and number of dimensions
        @testset "LatticeRule(z, s)" begin
            lattice_rule = LatticeRule(K_3600_32, 16)
            @test ndims(lattice_rule) == 16
            @test length(lattice_rule) == 2^32
        end

        # test constructor with generating vector only
        @testset "LatticeRule(z)" begin
            lattice_rule = LatticeRule(K_3600_32, 250)
            @test ndims(lattice_rule) == 250
            @test length(lattice_rule) == 2^32
        end

        # test constructor with number of dimensions only
        @testset "LatticeRule(s)" begin
            lattice_rule = LatticeRule(251)
            @test ndims(lattice_rule) == 251
            @test length(lattice_rule) == 2^32
        end

        # test constructor with file containing generating vector, number of dimensions and max number of points
        @testset "LatticeRule(z_file, s, n)" begin
            lattice_rule = LatticeRule(CKN_250_20_file, 9, 2^20)
            @test ndims(lattice_rule) == 9
            @test length(lattice_rule) == 2^20
        end

        # test constructor with file containing generating vector and number of dimensions
        @testset "LatticeRule(z_file, s)" begin
            lattice_rule = LatticeRule(K_3600_32_file, 16)
            @test ndims(lattice_rule) == 16
            @test length(lattice_rule) == 2^32
        end

        # test constructor with file containing generating vector only
        @testset "LatticeRule(z_file)" begin
            lattice_rule = LatticeRule(K_3600_32_file)
            @test ndims(lattice_rule) == 3600
            @test length(lattice_rule) == 2^32
        end

        # test getpoint
        @testset "getpoint(lattice_rule, k)" begin
            lattice_rule = LatticeRule(10)
            point0 = getpoint(lattice_rule, 0)
            @inferred getpoint(lattice_rule, 10)
            @test point0 isa eltype(LatticeRule)
            @test sum(point0) == 0
        end

        # test iterator access
        @testset "iterate(lattice_rule, state)" begin
            lattice_rule = LatticeRule(32)
            for i in lattice_rule
                nothing
            end
            @test all(first(lattice_rule) .== lattice_rule[0])
        end

        # test lattice_rule[i] access
        @testset "lattice_rule[i]" begin
            lattice_rule = LatticeRule(2)
            @test length(lattice_rule[0]) == 2
            @inferred lattice_rule[100]
            @test length(lattice_rule[end]) == 2
            @test firstindex(lattice_rule) == 0
            @test length(lattice_rule[1:20]) == 20
            @test length(collect(lattice_rule)) == 2^20
        end

        # test error handling
        @testset "error handling" begin
            lattice_rule = LatticeRule(9)
            @test_throws BoundsError getpoint(lattice_rule, -1)
            @test_throws BoundsError getpoint(lattice_rule, 2^20)
            @test_throws ArgumentError LatticeRule(0)
            @test_throws ArgumentError LatticeRule(-1)
            @test_throws ArgumentError LatticeRule(3601)
            @test_throws ArgumentError LatticeRule(CKN_250_20, 251)
            @test_throws ArgumentError LatticeRule(CKN_250_20_file, 12 ,0)
            @test_throws ArgumentError LatticeRule(K_3600_32, 251, 2^32 + 1)
        end

        # test print method
        @testset "show(lattice_rule)" begin
            lattice_rule = LatticeRule(9)
            str = string(lattice_rule)
        end

    end

    @testset "Constructing ShiftedLatticeRule" begin

        # test constructor with lattice rule and random shift
        @testset "ShiftedLatticeRule(lattice_rule, Δ)" begin
            lattice_rule = LatticeRule([UInt32(1), UInt32(5)], 2, 8)
            shifted_lattice_rule = ShiftedLatticeRule(lattice_rule, rand(2))
            @test ndims(shifted_lattice_rule) == 2
            @test length(shifted_lattice_rule) == 8
        end

        # test constructor with lattice rule only
        @testset "ShiftedLatticeRule(lattice_rule)" begin
            lattice_rule = LatticeRule([UInt32(1), UInt32(5)], 2, 8)
            shifted_lattice_rule = ShiftedLatticeRule(lattice_rule)
            @test ndims(shifted_lattice_rule) == 2
            @test length(shifted_lattice_rule) == 8
        end

        # test constructor with number of dimensions only
        @testset "ShiftedLatticeRule(s)" begin
            shifted_lattice_rule = ShiftedLatticeRule(251)
            @test ndims(shifted_lattice_rule) == 251
            @test length(shifted_lattice_rule) == 2^32
        end

        # test getpoint
        @testset "getpoint(shifted_lattice_rule, k)" begin
            shifted_lattice_rule = ShiftedLatticeRule(10)
            point0 = getpoint(shifted_lattice_rule, 0)
            @test point0 isa eltype(ShiftedLatticeRule)
            @inferred getpoint(shifted_lattice_rule, 101)
            @test all(0 .≤ point0 .≤ 1)
        end

        # test error handling
        @testset "error handling" begin
            lattice_rule = LatticeRule(100)
            @test_throws DimensionMismatch ShiftedLatticeRule(lattice_rule, rand(101))
            @test_throws ArgumentError ShiftedLatticeRule(lattice_rule, rand(100) .- 1)
            @test_throws ArgumentError ShiftedLatticeRule(lattice_rule, rand(100) .+ 1)
        end

        # test print method
        @testset "show(shifted_lattice_rule)" begin
            shifted_lattice_rule = ShiftedLatticeRule(9)
            str = string(shifted_lattice_rule)
        end

    end

    # approximate pi by throwing random darts in a square
    @testset "Approximating pi by throwing darts" begin
        darts(x) = x[1]*x[1] + x[2]*x[2] < 1
        lattice_rule = LatticeRule(2)
        Q = 4 * mean(darts.(collect(lattice_rule)))
        @test Q ≈ π rtol=1e-5
    end

    # see Keister, Bradley D. "Multidimensional quadrature algorithms." Computers in Physics 10.2 (1996): 119-128.
    @testset "Computing multidimensional integral from Keister, Bradley" begin
        dims = [9, 25, 60, 80, 100]
        exact = [-71.633234291 -1.356914e6 4.89052986e14 6.78878724e19 4.57024396e24]
        f(x) = cos(sqrt(sum(erfinv.(2*x .- 1).^2)))
        for (d, I) in Iterators.zip(dims, exact)
            lattice_rule = ShiftedLatticeRule(d)
            Q = π^(d/2) * mean(f.(collect(lattice_rule)))
            @test Q ≈ I rtol=1e-3
        end
    end

end
