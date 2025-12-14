# HyperbolicNumbers.jl

A Julia package for **hyperbolic (split-complex, perplex or dual) numbers** of the form

[ z = a + b hy, hy^2 = 1 with hy ≠ ±1]

This package provides a mathematically rigorous and Julia-idiomatic implementation of hyperbolic numbers, including algebraic operations, geometric (Minkowski) interpretation, and plotting support.

---

## ✨ Features

* Fully implemented `Number` interface
* Addition, subtraction, multiplication, division, and inversion
* Hyperbolic conjugation and quadratic norm
* Correct handling of **zero divisors**
* Minkowski (Lorentzian) geometry interpretation
* Classification into timelike, spacelike, and lightlike elements
* Hyperbolic angle (rapidity)
* Plot recipes compatible with `Plots.jl`
* Extensive inline documentation and docstrings

---

## 📦 Installation

From the Julia REPL:

```julia
pkg> add https://github.com/elletri/HyperbolicNumbers.jl
```

For local development:

```julia
pkg> dev path="/path/to/HyperbolicNumbers"
```

---

## 🚀 Quick start

```julia
using HyperbolicNumbers

z1 = Hyperbolic(2, 1)
z2 = Hyperbolic(1, -1)

z1 + z2
z1 * z2
inv(z1)
conj(z1)
abs(z1)  # a^2 - b^2
```

Mixed arithmetic with real numbers works naturally:

```julia
z = Hyperbolic(3, 2)
2z + 1
```

---

## 🧮 Algebraic properties

### Hyperbolic unit

```julia
hy^2 == one(Hyperbolic{Int})
```

### Zero divisors

Hyperbolic numbers are **not a field**. Elements with

[ a^2 - b^2 = 0 ]

are nonzero but non-invertible.

```julia
z = Hyperbolic(1, 1)
islightlike(z)   # true
isunit(z)        # false
```

---

## 🌌 Minkowski geometry

Each hyperbolic number corresponds to a vector in 2D Minkowski space.

```julia
z = Hyperbolic(3, 1)

istimelike(z)
isspacelike(z)
islightlike(z)

lorentz_inner(z, z)
hyperbolic_angle(z)
```

The quadratic form

[ \lVert z \rVert = a^2 - b^2 ]

is the Lorentzian norm.

---

## 📈 Visualization

Hyperbolic numbers can be plotted directly using `Plots.jl`:

```julia
using Plots, HyperbolicNumbers

z = Hyperbolic(2, 1)
plot(z)
```

This plots the point `(a, b)` in the Minkowski plane.

For low-level control:

```julia
x, y = visualize(z)
scatter([x], [y])
```

---

## 📚 Mathematical background

Hyperbolic numbers (also called *split-complex*, *perplex*, or *Lorentz numbers*) form a two-dimensional commutative real algebra with basis `{1, h}` and defining relation `h² = 1`.

They naturally encode:

* 2D special relativity
* Lorentz transformations
* Hyperbolic rotations
* Indefinite inner product spaces

---

## 🧪 Testing

Run tests with:

```julia
pkg> test HyperbolicNumbers
```

---

## 🤝 Contributing

Contributions are welcome. Please:

* add tests for new features
* document all public functions
* keep mathematical behavior explicit

---

## 🙏 Acknowledgements

This package was developed with the assistance of OpenAI's ChatGPT, which was used as a programming and documentation aid during development. All design choices, mathematical interpretations, and any remaining errors are the sole responsibility of the author.

Info on hyperbolic numbers in Julia is also in https://doi.org/10.48550/arXiv.2301.01707

---

## 📄 License

MIT License

---

## 🔗 Related work

* Complex numbers (`Complex{T}` in Base Julia)
* Dual numbers
* Clifford algebras
* Lorentzian geometry

---

**Author:** Lorenzo Tripodi
**Status:** Experimental / Research-ready
