using Test
using HyperbolicNumbers
using Plots

@testset "HyperbolicNumbers" begin
    

    #println("=== HyperbolicNumbers.jl test suite ===")

    # --------------------------------------------------
    # Construction
    # --------------------------------------------------
    z = Hyperbolic(3.0, 2.0)
    w = Hyperbolic(2.0, -1.0)   # invertible: 4 - 1 ≠ 0

    @test z.a == 3.0
    @test z.b == 2.0

    #println("✔ Construction")

    # --------------------------------------------------
    # Basic arithmetic
    # --------------------------------------------------
    @test z + w == Hyperbolic(5.0, 1.0)
    @test z - w == Hyperbolic(1.0, 3.0)
    @test z * w == Hyperbolic(4.0, 1.0)

    #println("✔ Arithmetic")

    # --------------------------------------------------
    # Division (invertible only)
    # --------------------------------------------------
    @test z / w ≈ z * inv(w)
    @test w / w ≈ one(w)

    #println("✔ Division")

    # --------------------------------------------------
    # Mixed arithmetic with reals
    # --------------------------------------------------
    @test z + 1 == Hyperbolic(4.0, 2.0)
    @test 2 * z == Hyperbolic(6.0, 4.0)
    @test z - 2 == Hyperbolic(1.0, 2.0)

    #println("✔ Mixed arithmetic")

    # --------------------------------------------------
    # Conjugation and norm
    # --------------------------------------------------
    zc = conj(z)
    @test zc == Hyperbolic(3.0, -2.0)
    @test abs(z) == 5.0

    #println("✔ Conjugation and quadratic norm")

    # --------------------------------------------------
    # Inversion and units
    # --------------------------------------------------
    @test isunit(z)
    @test inv(z) * z ≈ one(z)

    #println("✔ Inversion and units")

    # --------------------------------------------------
    # Zero divisors
    # --------------------------------------------------
    l = Hyperbolic(1.0, 1.0)

    @test islightlike(l)
    @test !isunit(l)

    #println("✔ Zero divisors")

    # --------------------------------------------------
    # Minkowski classification
    # --------------------------------------------------
    t = Hyperbolic(3.0, 1.0)
    s = Hyperbolic(1.0, 3.0)

    @test istimelike(t)
    @test isspacelike(s)
    @test !islightlike(t)

    #println("✔ Minkowski classification")

    # --------------------------------------------------
    # Lorentz inner product
    # --------------------------------------------------
    @test lorentz_inner(z, z) == abs(z)

    #println("✔ Lorentz inner product")

    # --------------------------------------------------
    # Hyperbolic angle
    # --------------------------------------------------
    θ = hyperbolic_angle(t)
    @test isfinite(θ)

    #println("✔ Hyperbolic angle")

    # --------------------------------------------------
    # Visualization helper
    # --------------------------------------------------
    x, y = visualize(z)
    @test x == z.a
    @test y == z.b

    #println("✔ Visualization helper")

    # --------------------------------------------------
    # Plot recipe (smoke test)
    # --------------------------------------------------
    #p = plot(z)
    @test_nowarn p = plot(z)
    @test p isa Plots.Plot
    
    #println("✔ Plot recipe")

    # --------------------------------------------------
    # Expected failures
    # --------------------------------------------------
    try
        inv(l)
        error("inv(l) should have failed for zero divisor")
    catch e
        @test e isa DomainError
    end

    #println("✔ Expected failures")

    println("=== All HyperbolicNumbers tests passed ===")
end
