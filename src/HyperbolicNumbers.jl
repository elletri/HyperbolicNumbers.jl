module HyperbolicNumbers

# --------------------------------------
# Julia package to support hyperbolic numbers in Julia
# Created with the support of ChatGPT
#
# L. Tripodi, EHV, Sat 13 Dec 2025
# ---------------------------------------



export Hyperbolic, h, conj, abs, visualize,
       inv, isunit, islightlike, istimelike, isspacelike,
       lorentz_inner, hyperbolic_angle

using RecipesBase

"""
    Hyperbolic{T<:Real}(a, b)

Represents a **hyperbolic (split-complex) number** of the form

```
a + h*b
```

where `a` and `b` are real numbers and `h` is the *hyperbolic unit* satisfying

```
h^2 = 1
```

Hyperbolic numbers form a commutative algebra over the reals and implement
Julia's `Number` interface. They naturally encode **2D Minkowski geometry**.

# Fields
- `a::T`: real (scalar) part
- `b::T`: hyperbolic part
"""
struct Hyperbolic{T<:Real} <: Number
    a::T
    b::T
end

"""
    h

The **hyperbolic unit**, defined such that `h^2 = 1`.
"""
const h = Hyperbolic(0, 1)

"""
    Hyperbolic(a::Real)

Constructs a purely real hyperbolic number `a + h*0`.
"""
Hyperbolic(a::Real) = Hyperbolic(a, zero(a))

# ---------------------------------------------------------------------
# Number interface
# ---------------------------------------------------------------------

"""
    zero(::Type{Hyperbolic{T}})

Additive identity `0 + h*0`.
"""
Base.zero(::Type{Hyperbolic{T}}) where T = Hyperbolic(zero(T), zero(T))

"""
    one(::Type{Hyperbolic{T}})

Multiplicative identity `1 + h*0`.
"""
Base.one(::Type{Hyperbolic{T}}) where T = Hyperbolic(one(T), zero(T))

"""
    real(z::Hyperbolic)

Returns the real part `a`.
"""
Base.real(z::Hyperbolic) = z.a

"""
    imag(z::Hyperbolic)

Returns the hyperbolic component `b`.
"""
Base.imag(z::Hyperbolic) = z.b

"""
    iszero(z::Hyperbolic)

Checks whether `z == 0 + h*0`.
"""
Base.iszero(z::Hyperbolic) = iszero(z.a) && iszero(z.b)

"""
    promote_rule(::Type{Hyperbolic{T}}, ::Type{S})

Allows mixed arithmetic with real numbers.
"""
Base.promote_rule(::Type{Hyperbolic{T}}, ::Type{S}) where {T<:Real,S<:Real} = Hyperbolic{promote_type(T,S)}

# ---------------------------------------------------------------------
# Display and comparison
# ---------------------------------------------------------------------

"""
    show(io, z::Hyperbolic)

Pretty-print as `a + h*b`.
"""
Base.show(io::IO, z::Hyperbolic) = print(io, "$(z.a) + h*$(z.b)")

"""
    ==(z1::Hyperbolic, z2::Hyperbolic)

Component-wise equality.
"""
Base.:(==)(z1::Hyperbolic, z2::Hyperbolic) = z1.a == z2.a && z1.b == z2.b

# ---------------------------------------------------------------------
# Arithmetic
# ---------------------------------------------------------------------

"""
    +(z1::Hyperbolic, z2::Hyperbolic)

Component-wise addition.
"""
Base.:+(z1::Hyperbolic, z2::Hyperbolic) = Hyperbolic(z1.a + z2.a, z1.b + z2.b)

"""
    -(z1::Hyperbolic, z2::Hyperbolic)

Component-wise subtraction.
"""
Base.:-(z1::Hyperbolic, z2::Hyperbolic) = Hyperbolic(z1.a - z2.a, z1.b - z2.b)

"""
    -(z::Hyperbolic)

Unary negation.
"""
Base.:-(z::Hyperbolic) = Hyperbolic(-z.a, -z.b)

"""
    *(z1::Hyperbolic, z2::Hyperbolic)

Hyperbolic multiplication using `h^2 = 1`.
"""
Base.:*(z1::Hyperbolic, z2::Hyperbolic) =
    Hyperbolic(z1.a*z2.a + z1.b*z2.b,
               z1.a*z2.b + z1.b*z2.a)

"""
    conj(z::Hyperbolic)

Hyperbolic conjugation: `a - h*b`.
"""
conj(z::Hyperbolic) = Hyperbolic(z.a, -z.b)

"""
    abs(z::Hyperbolic)

Quadratic (Minkowski) norm: `a^2 - b^2`.
"""
abs(z::Hyperbolic) = z.a^2 - z.b^2

"""
    inv(z::Hyperbolic)

Multiplicative inverse of `z`.

Throws `DivideError` if `abs(z) == 0` (zero divisors).
"""
Base.inv(z::Hyperbolic) = begin
    n = abs(z)
    n == 0 && throw(DivideError())
    w = conj(z)
    Hyperbolic(w.a / n, w.b / n)
end

"""
    /(z1::Hyperbolic, z2::Hyperbolic)

Division via `z1 * inv(z2)`.
"""
Base.:/(z1::Hyperbolic, z2::Hyperbolic) = z1 * inv(z2)

"""
    isunit(z::Hyperbolic)

Returns `true` if `z` is invertible, i.e. `abs(z) ≠ 0`.
"""
isunit(z::Hyperbolic) = abs(z) != 0

# ---------------------------------------------------------------------
# Geometry and classification
# ---------------------------------------------------------------------

"""
    islightlike(z::Hyperbolic)

Returns `true` if `abs(z) == 0` and `z ≠ 0`.
"""
islightlike(z::Hyperbolic) = abs(z) == 0 && !iszero(z)

"""
    istimelike(z::Hyperbolic)

Returns `true` if `abs(z) > 0`.
"""
istimelike(z::Hyperbolic) = abs(z) > 0

"""
    isspacelike(z::Hyperbolic)

Returns `true` if `abs(z) < 0`.
"""
isspacelike(z::Hyperbolic) = abs(z) < 0

"""
    lorentz_inner(z1::Hyperbolic, z2::Hyperbolic)

Lorentzian inner product:

```
⟨z1, z2⟩ = a₁a₂ − b₁b₂
```
"""
lorentz_inner(z1::Hyperbolic, z2::Hyperbolic) = z1.a*z2.a - z1.b*z2.b

"""
    hyperbolic_angle(z::Hyperbolic)

Returns the hyperbolic angle (rapidity) θ such that

```
z = r * (cosh(θ) + h*sinh(θ))
```

Defined only for timelike elements.
"""
hyperbolic_angle(z::Hyperbolic) = begin
    istimelike(z) || throw(DomainError(z, "Hyperbolic angle defined only for timelike elements"))
    atanh(z.b / z.a)
end

# ---------------------------------------------------------------------
# Elementary functions
# ---------------------------------------------------------------------

"""
    exp(z::Hyperbolic)

Exponential of a hyperbolic number.
"""
Base.exp(z::Hyperbolic) = begin
    ea = exp(z.a)
    Hyperbolic(ea*cosh(z.b), ea*sinh(z.b))
end

"""
    sinh(z::Hyperbolic)

Hyperbolic sine.
"""
Base.sinh(z::Hyperbolic) = Hyperbolic(sinh(z.a)*cosh(z.b), cosh(z.a)*sinh(z.b))

"""
    cosh(z::Hyperbolic)

Hyperbolic cosine.
"""
Base.cosh(z::Hyperbolic) = Hyperbolic(cosh(z.a)*cosh(z.b), sinh(z.a)*sinh(z.b))

# ---------------------------------------------------------------------
# Visualization and plotting
# ---------------------------------------------------------------------

"""
    visualize(z::Hyperbolic)

Returns `(a, b)` for plotting in the Minkowski plane.
"""
visualize(z::Hyperbolic) = (z.a, z.b)

"""
Plot recipe for hyperbolic numbers.

Allows direct plotting via:

```julia
plot(Hyperbolic(2,1))
```
"""
@recipe function f(z::Hyperbolic)
    seriestype := :scatter
    xlabel --> "a (real part)"
    ylabel --> "b (hyperbolic part)"
    [z.a], [z.b]
end

end # module

