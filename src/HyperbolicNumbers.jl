module HyperbolicNumbers

# --------------------------------------
# Julia package to support hyperbolic numbers (also known as perplex numbers, split-complex numbers or double numbers)
# Created with the support of ChatGPT
#
# L. Tripodi, EHV, Sat 13 Dec 2025
# ---------------------------------------

# Changes
# Added alternative names of the hyperbolic numbers in the title
# Corrected the visualization function to resemble more the visualization of complex numbers in Julia. 
# The hyperbolic unit is hy a bit as the imaginary unit is im in Julia

export Hyperbolic, hy, conj, abs, visualize,
       inv, isunit, islightlike, istimelike, isspacelike,
       lorentz_inner, hyperbolic_angle

using RecipesBase

"""
    Hyperbolic{T<:Real}(a, b)

Represents a **hyperbolic (split-complex, perplex, double) number** of the form

```
a + b*hy
```

where `a` and `b` are real numbers and `hy` is the *hyperbolic unit* satisfying

```
hy^2 = 1
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
    Hyperbolic(a::Real, b::Real)

Construct a hyperbolic number `a + h*b` with automatic promotion
of `a` and `b` to the same type.
"""
Hyperbolic(a::Real, b::Real) = Hyperbolic(promote(a, b)...)


"""
    hy

The **hyperbolic unit**, defined such that `hy^2 = 1`.
"""
const hy = Hyperbolic(0, 1)

"""
    Hyperbolic(a::Real)

Constructs a purely real hyperbolic number `a + 0*hy`.
"""
Hyperbolic(a::Real) = Hyperbolic(a, zero(a))

# ---------------------------------------------------------------------
# Number interface
# ---------------------------------------------------------------------

"""
    zero(::Type{Hyperbolic{T}})

Additive identity `0 + hy*0`.
"""
Base.zero(::Type{Hyperbolic{T}}) where T = Hyperbolic(zero(T), zero(T))

"""
    one(::Type{Hyperbolic{T}})

Multiplicative identity `1 + hy*0`.
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
Base.promote_rule(::Type{S}, ::Type{Hyperbolic{T}}) where {T<:Real,S<:Real} = Hyperbolic{promote_type(T,S)}

# ---------------------------------------------------------------------
# Display and comparison
# ---------------------------------------------------------------------

"""
    show(io, z::Hyperbolic)

Pretty-print as `a + b*hy`.
"""
Base.show(io::IO, z::Hyperbolic) = print(io, "$(z.a) + $(z.b)hy")

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

Hyperbolic multiplication using `hy^2 = 1`.
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

Return the multiplicative inverse of `z`.
Throws `DomainError` if `z` is a zero divisor.
"""
function Base.inv(z::Hyperbolic{T}) where {T}
    n = abs(z)
    if iszero(n)
        throw(DomainError(z, "Hyperbolic number is not invertible (zero divisor)"))
    end
    Hyperbolic(z.a/n, -z.b/n)
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

# -----------------------------
# Mixed arithmetic with reals
# -----------------------------

Base.:+(z::Hyperbolic, x::Real) = z + Hyperbolic(x)
Base.:+(x::Real, z::Hyperbolic) = Hyperbolic(x) + z
Base.:-(z::Hyperbolic, x::Real) = z - Hyperbolic(x)
Base.:-(x::Real, z::Hyperbolic) = Hyperbolic(x) - z
Base.:*(z::Hyperbolic, x::Real) = Hyperbolic(z.a*x, z.b*x)
Base.:*(x::Real, z::Hyperbolic) = z * x
Base.:/(z::Hyperbolic, x::Real) = Hyperbolic(z.a/x, z.b/x)


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

