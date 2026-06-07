# chern-classes

> **Characteristic classes for vector bundles — the bridge between geometry and topology.**

[![crates.io](https://img.shields.io/crates/v/chern-classes.svg)](https://crates.io/crates/chern-classes)
[![docs.rs](https://docs.rs/chern-classes/badge.svg)](https://docs.rs/chern-classes)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![tests](https://img.shields.io/badge/tests-16-passing-green.svg)]()

Computes Chern classes, Pontryagin classes, Todd classes, and related invariants for complex and real vector bundles. Implements the splitting principle for decomposing bundles into line bundles.

---

## Why This Exists

Characteristic classes are the **central tool** of differential topology — they translate geometric information (curvature, connections) into topological invariants (cohomology classes). But working with them typically requires:

- Dense textbooks (Milnor & Stasheff, 300+ pages)
- Symbolic algebra systems (Mathematica, Sage)
- Manual polynomial manipulation

`chern-classes` makes characteristic classes **computable** — a Rust library that handles the polynomial algebra so you can focus on the geometry.

---

## Architecture

```
                    ┌─────────────────────┐
                    │   Vector Bundle E    │
                    │   rank n, type C/R   │
                    └────────┬────────────┘
                             │
              ┌──────────────┼──────────────┐
              │              │              │
     ┌────────▼───┐  ┌──────▼─────┐  ┌─────▼──────────┐
     │ Chern      │  │ Pontryagin │  │ Todd           │
     │ Classes    │  │ Classes    │  │ Class          │
     │ c₁...cₙ   │  │ p₁...pₙ/₂ │  │ td(E)         │
     └────────┬───┘  └──────┬─────┘  └─────┬──────────┘
              │              │              │
              └──────┬───────┘              │
                     │                      │
          ┌──────────▼──────────┐           │
          │ Splitting Principle │◄──────────┘
          │ E → L₁ ⊕ ... ⊕ Lₙ │
          └─────────────────────┘
                     │
          ┌──────────▼──────────┐
          │ Polynomial Ring     │
          │ (add, mul, eval)    │
          └─────────────────────┘
```

---

## Installation

```toml
[dependencies]
chern-classes = "0.1.0"
```

---

## Quick Start

```rust
use chern_classes::{ChernCalculator, Polynomial};

// Chern class of a line bundle with c₁ = 3
let c = ChernCalculator::line_bundle_class(3.0);
// c = 1 + 3x

// Whitney sum: c(E ⊕ F) = c(E) · c(F)
let c1 = ChernCalculator::line_bundle_class(1.0);
let c2 = ChernCalculator::line_bundle_class(2.0);
let total = ChernCalculator::whitney_sum(&c1, &c2);
// (1+x)(1+2x) = 1 + 3x + 2x²
assert!((total.coeff(1) - 3.0).abs() < 0.001);
assert!((total.coeff(2) - 2.0).abs() < 0.001);
```

---

## Usage Examples

### Example 1: Chern Character

```rust
use chern_classes::ChernCalculator;

// Chern character: ch(E) = rank + c₁ + (c₁² - 2c₂)/2
let ch = ChernCalculator::chern_character(2, 3.0, 1.0);
assert!((ch.coeff(0) - 2.0).abs() < 0.001);  // rank
assert!((ch.coeff(1) - 3.0).abs() < 0.001);  // c₁
assert!((ch.coeff(2) - 3.5).abs() < 0.001);  // (9-2)/2 = 3.5
```

### Example 2: Todd Class (for Riemann-Roch)

```rust
use chern_classes::ToddClass;

// td(E) = 1 + c₁/2 + (c₁² + c₂)/12
let td = ToddClass::compute(2.0, 1.0);
assert!((td.coeff(0) - 1.0).abs() < 0.001);     // 1
assert!((td.coeff(1) - 1.0).abs() < 0.001);     // c₁/2 = 1
assert!((td.coeff(2) - 5.0/12.0).abs() < 0.01); // (4+1)/12
```

### Example 3: Splitting Principle

```rust
use chern_classes::SplittingPrinciple;

// Decompose rank-2 bundle with Chern roots x₁=1, x₂=2
let total = SplittingPrinciple::split(2, &[1.0, 2.0]);
// (1+x)(1+2x) = 1 + 3x + 2x²
assert!((total.coeff(0) - 1.0).abs() < 0.001);
assert!((total.coeff(1) - 3.0).abs() < 0.001);

// Elementary symmetric polynomials
let sigma = SplittingPrinciple::elementary_symmetric(&[1.0, 2.0, 3.0]);
// σ₁ = 6, σ₂ = 11, σ₃ = 6
assert!((sigma[1] - 6.0).abs() < 0.001);
assert!((sigma[2] - 11.0).abs() < 0.001);
assert!((sigma[3] - 6.0).abs() < 0.001);
```

### Example 4: Pontryagin Classes

```rust
use chern_classes::PontryaginClass;

// p₁(E) = c₁² for real bundles
let p1 = PontryaginClass::first(2.0);
assert!((p1.coeff(2) - 4.0).abs() < 0.001);

// Total Pontryagin class
let p = PontryaginClass::total(1.0, 0.5);
assert!((p.coeff(0) - 1.0).abs() < 0.001);
assert!((p.coeff(2) - 1.0).abs() < 0.001);
```

---

## Mathematical Background

### Chern Classes

For a complex vector bundle E → M of rank n:

- **c₀(E) = 1** (by convention)
- **c₁(E)** ∈ H²(M) — first Chern class
- **cₖ(E)** ∈ H²ᵏ(M) — k-th Chern class
- **c(E) = 1 + c₁ + ... + cₙ** — total Chern class

**Whitney sum formula**: c(E ⊕ F) = c(E) · c(F)

### Todd Class

The Todd class appears in the **Hirzebruch-Riemann-Roch theorem**:

χ(E) = ∫ ch(E) · td(TM)

For a line bundle L with first Chern class x:
td(L) = x/(1 - e⁻ˣ) = 1 + x/2 + x²/12 + ...

### Pontryagin Classes

For a real vector bundle E, Pontryagin classes are:
pₖ(E) = (-1)ᵏ c₂ₖ(E ⊗ ℂ) ∈ H⁴ᵏ(M)

### Splitting Principle

The splitting principle states that any rank-n bundle can be formally decomposed as a direct sum of line bundles: E = L₁ ⊕ ... ⊕ Lₙ, with Chern classes expressed in terms of the Chern roots x₁, ..., xₙ.

---

## API Reference

| Type | Description |
|------|-------------|
| `Polynomial` | Polynomial ring with add, mul, evaluate, display |
| `VectorBundle` | Named bundle with rank and properties |
| `ChernCalculator` | Chern classes, Whitney sum, Chern character |
| `ToddClass` | Todd class computation |
| `PontryaginClass` | Real bundle characteristic classes |
| `SplittingPrinciple` | Decompose bundles into line bundles |

---

## Performance

- **Zero dependencies** — no external math libraries
- **F64 precision** — suitable for numerical computations
- **Stack-allocated polynomials** — no heap allocation for typical use

```
test result: ok. 16 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
```

---

## References

- Milnor, J. & Stasheff, J. *Characteristic Classes* (1974)
- Hirzebruch, F. *Topological Methods in Algebraic Geometry* (1966)
- Bott, R. & Tu, L. *Differential Forms in Algebraic Topology* (1982)

---

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project.*
