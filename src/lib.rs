//! # Chern Classes
//!
//! A library for computing characteristic classes of vector bundles.
//!
//! Provides:
//! - Vector bundle representation
//! - Chern classes with Whitney sum formula
//! - Todd class for Riemann-Roch theorem
//! - Pontryagin classes for real vector bundles
//! - Splitting principle for decomposing bundles

/// A polynomial in one variable with rational coefficients.
/// Stored as a vector of coefficients: coeffs[i] is the coefficient of t^i.
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    coeffs: Vec<f64>,
}

impl Polynomial {
    /// Create a polynomial from coefficients.
    pub fn new(coeffs: Vec<f64>) -> Self {
        let mut p = Self { coeffs };
        p.trim();
        p
    }

    /// The zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: vec![] }
    }

    /// The constant polynomial 1.
    pub fn one() -> Self {
        Self { coeffs: vec![1.0] }
    }

    /// Create a monomial c·t^n.
    pub fn monomial(coeff: f64, degree: usize) -> Self {
        if coeff == 0.0 { return Self::zero(); }
        let mut coeffs = vec![0.0; degree + 1];
        coeffs[degree] = coeff;
        Self { coeffs }
    }

    /// Degree of the polynomial.
    pub fn degree(&self) -> usize {
        if self.coeffs.is_empty() { 0 } else { self.coeffs.len() - 1 }
    }

    /// Evaluate at a point.
    pub fn evaluate(&self, t: f64) -> f64 {
        let mut result = 0.0;
        let mut power = 1.0;
        for &c in &self.coeffs {
            result += c * power;
            power *= t;
        }
        result
    }

    /// Add two polynomials.
    pub fn add(&self, other: &Polynomial) -> Polynomial {
        let len = self.coeffs.len().max(other.coeffs.len());
        let mut coeffs = vec![0.0; len];
        for (i, &c) in self.coeffs.iter().enumerate() {
            coeffs[i] += c;
        }
        for (i, &c) in other.coeffs.iter().enumerate() {
            coeffs[i] += c;
        }
        Polynomial::new(coeffs)
    }

    /// Subtract polynomials.
    pub fn sub(&self, other: &Polynomial) -> Polynomial {
        let neg = Polynomial::new(other.coeffs.iter().map(|&c| -c).collect());
        self.add(&neg)
    }

    /// Multiply two polynomials.
    pub fn mul(&self, other: &Polynomial) -> Polynomial {
        if self.coeffs.is_empty() || other.coeffs.is_empty() {
            return Polynomial::zero();
        }
        let len = self.coeffs.len() + other.coeffs.len() - 1;
        let mut coeffs = vec![0.0; len];
        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in other.coeffs.iter().enumerate() {
                coeffs[i + j] += a * b;
            }
        }
        Polynomial::new(coeffs)
    }

    /// Scale by a constant.
    pub fn scale(&self, scalar: f64) -> Polynomial {
        Polynomial::new(self.coeffs.iter().map(|&c| c * scalar).collect())
    }

    /// Get the coefficient of degree k.
    pub fn coeff(&self, k: usize) -> f64 {
        self.coeffs.get(k).copied().unwrap_or(0.0)
    }

    /// Formal derivative.
    pub fn derivative(&self) -> Polynomial {
        if self.degree() == 0 { return Polynomial::zero(); }
        let coeffs: Vec<f64> = self.coeffs.iter().skip(1)
            .enumerate()
            .map(|(i, &c)| c * (i + 1) as f64)
            .collect();
        Polynomial::new(coeffs)
    }

    fn trim(&mut self) {
        while self.coeffs.last() == Some(&0.0) {
            self.coeffs.pop();
        }
    }
}

impl std::fmt::Display for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coeffs.is_empty() {
            return write!(f, "0");
        }
        let mut first = true;
        for (i, &c) in self.coeffs.iter().enumerate() {
            if c == 0.0 { continue; }
            if !first { write!(f, " + ")?; }
            if i == 0 {
                write!(f, "{}", c)?;
            } else if i == 1 {
                write!(f, "{}t", c)?;
            } else {
                write!(f, "{}t^{}", c, i)?;
            }
            first = false;
        }
        if first { write!(f, "0")?; }
        Ok(())
    }
}

/// A vector bundle.
#[derive(Debug, Clone)]
pub struct VectorBundle {
    /// Rank of the bundle.
    pub rank: usize,
    /// Name/label.
    pub name: String,
    /// Chern classes c_0, c_1, ..., c_rank as polynomials in the base cohomology.
    /// c_0 = 1 always.
    chern_classes: Vec<Polynomial>,
}

impl VectorBundle {
    /// Create a trivial bundle of given rank.
    pub fn trivial(rank: usize) -> Self {
        let chern = vec![Polynomial::one()];
        // Trivial bundle has c_k = 0 for k > 0
        Self { rank, name: format!("O^{}", rank), chern_classes: chern }
    }

    /// Create a line bundle (rank 1) with given first Chern class.
    pub fn line_bundle(c1: Polynomial) -> Self {
        Self {
            rank: 1,
            name: "L".to_string(),
            chern_classes: vec![Polynomial::one(), c1],
        }
    }

    /// Create the tautological line bundle O(-1) over CP^n.
    pub fn tautological_line_bundle() -> Self {
        // c_1(O(-1)) = -h where h is the hyperplane class
        Self::line_bundle(Polynomial::monomial(-1.0, 1))
    }

    /// Create the hyperplane line bundle O(1) over CP^n.
    pub fn hyperplane_bundle() -> Self {
        Self::line_bundle(Polynomial::monomial(1.0, 1))
    }

    /// Create the tangent bundle of CP^n.
    pub fn tangent_bundle_cp(n: usize) -> Self {
        // T(CP^n) ⊕ O = O(1)^{⊕(n+1)}
        // c(T) = (1 + h)^{n+1}
        let h = Polynomial::new(vec![1.0, 1.0]); // 1 + h
        let mut total = Polynomial::one();
        for _ in 0..=n {
            total = total.mul(&h);
        }
        // Expand into individual Chern classes
        let mut chern = Vec::new();
        for k in 0..=n {
            chern.push(Polynomial::monomial(total.coeff(k), k));
        }
        Self {
            rank: n,
            name: format!("T(CP^{})", n),
            chern_classes: chern,
        }
    }

    /// Get the k-th Chern class as a reference (panics if out of bounds for k=0).
    /// Use get_chern_class for safe access.
    pub fn chern_class(&self, k: usize) -> &Polynomial {
        &self.chern_classes[k]
    }

    /// Get the total Chern class c = 1 + c_1 + c_2 + ...
    pub fn total_chern_class(&self) -> Polynomial {
        let mut total = Polynomial::zero();
        for (_k, ck) in self.chern_classes.iter().enumerate() {
            total = total.add(ck);
        }
        total
    }

    /// Get all Chern classes.
    pub fn chern_classes(&self) -> &[Polynomial] {
        &self.chern_classes
    }
}

impl VectorBundle {
    /// Get the k-th Chern class (owned, safe).
    pub fn get_chern_class(&self, k: usize) -> Polynomial {
        if k < self.chern_classes.len() {
            self.chern_classes[k].clone()
        } else {
            Polynomial::zero()
        }
    }
}

/// Chern class computations.
pub struct ChernClass;

impl ChernClass {
    /// Whitney sum formula: c(E ⊕ F) = c(E) · c(F).
    pub fn whitney_sum(e: &VectorBundle, f: &VectorBundle) -> VectorBundle {
        let rank = e.rank + f.rank;
        let c_e = e.total_chern_class();
        let c_f = f.total_chern_class();
        let c_total = c_e.mul(&c_f);

        // Expand into individual classes
        let mut chern = vec![Polynomial::one()];
        for k in 1..=rank {
            chern.push(Polynomial::monomial(c_total.coeff(k), k));
        }

        VectorBundle {
            rank,
            name: format!("{} ⊕ {}", e.name, f.name),
            chern_classes: chern,
        }
    }

    /// Compute the top Chern class c_n of a rank-n bundle.
    pub fn top_class(bundle: &VectorBundle) -> Polynomial {
        bundle.get_chern_class(bundle.rank)
    }

    /// Compute c_1 for a tensor product of line bundles.
    /// c_1(L₁ ⊗ L₂) = c_1(L₁) + c_1(L₂).
    pub fn c1_tensor_line(l1: &VectorBundle, l2: &VectorBundle) -> Polynomial {
        assert_eq!(l1.rank, 1);
        assert_eq!(l2.rank, 1);
        l1.get_chern_class(1).add(&l2.get_chern_class(1))
    }

    /// Compute the Chern character ch(E) = rank + c_1 + (c_1² - 2c_2)/2 + ...
    pub fn chern_character(bundle: &VectorBundle) -> Polynomial {
        let rank = Polynomial::monomial(bundle.rank as f64, 0);
        let c1 = bundle.get_chern_class(1);
        let c2 = bundle.get_chern_class(2);

        // ch = rank + c1 + (c1² - 2c2)/2
        let c1_sq = c1.mul(&c1);
        let two_c2 = c2.scale(2.0);
        let ch2 = c1_sq.sub(&two_c2).scale(0.5);

        rank.add(&c1).add(&ch2)
    }
}

/// Todd class computation for the Hirzebruch-Riemann-Roch theorem.
pub struct ToddClass;

impl ToddClass {
    /// Compute the Todd class from Chern classes.
    /// For a line bundle with first Chern class x:
    ///   td = x / (1 - e^{-x})
    /// For small x, td ≈ 1 + x/2 + x²/12 + ...
    pub fn compute(bundle: &VectorBundle) -> Polynomial {
        let c1 = bundle.get_chern_class(1);
        let c2 = bundle.get_chern_class(2);
        let c1_sq = c1.mul(&c1);

        // td = 1 + c1/2 + (c1² + c2)/12 + c1·c2/24 + ...
        let td1 = c1.scale(0.5);
        let td2 = c1_sq.add(&c2).scale(1.0 / 12.0);
        let td3 = c1.mul(&c2).scale(1.0 / 24.0);

        Polynomial::one().add(&td1).add(&td2).add(&td3)
    }

    /// Compute the Todd genus (integral of the Todd class over the manifold).
    /// For CP^n, the Todd genus is 1.
    pub fn todd_genus(bundle: &VectorBundle) -> f64 {
        let td = Self::compute(bundle);
        // The top component of the Todd class
        td.coeff(bundle.rank)
    }
}

/// Pontryagin classes for real vector bundles.
///
/// For a real vector bundle E of rank n, the k-th Pontryagin class is
/// p_k(E) = (-1)^k c_{2k}(E ⊗ C) ∈ H^{4k}(M).
pub struct PontryaginClass;

impl PontryaginClass {
    /// Compute Pontryagin classes from Chern classes of the complexified bundle.
    pub fn from_chern(bundle: &VectorBundle) -> Vec<Polynomial> {
        let mut pontryagin = vec![Polynomial::one()];

        for k in 1..=(bundle.rank / 2) {
            let c2k = bundle.get_chern_class(2 * k);
            // p_k = (-1)^k c_{2k}
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            pontryagin.push(c2k.scale(sign));
        }
        pontryagin
    }

    /// The first Pontryagin class p_1.
    pub fn p1(bundle: &VectorBundle) -> Polynomial {
        let c2 = bundle.get_chern_class(2);
        c2.scale(-1.0) // p_1 = -c_2
    }

    /// The Â-genus from Pontryagin classes.
    /// Â = 1 - p_1/24 + (7p_1² - 4p_2)/5760 + ...
    pub fn a_genus(bundle: &VectorBundle) -> Polynomial {
        let p1 = Self::p1(bundle);
        let p1_sq = p1.mul(&p1);

        let term1 = p1.scale(-1.0 / 24.0);
        let term2 = p1_sq.scale(7.0 / 5760.0);

        Polynomial::one().add(&term1).add(&term2)
    }
}

/// Splitting principle.
///
/// The splitting principle allows us to treat a rank-n vector bundle E
/// as a direct sum of line bundles L₁ ⊕ ... ⊕ Lₙ after pulling back
/// to a suitable space.
pub struct SplittingPrinciple;

impl SplittingPrinciple {
    /// Decompose a rank-n bundle into formal Chern roots x₁, ..., xₙ.
    /// The Chern classes are the elementary symmetric polynomials in x₁, ..., xₙ.
    pub fn chern_roots(rank: usize) -> Vec<Polynomial> {
        (0..rank).map(|i| Polynomial::monomial(1.0, 1).add(&Polynomial::monomial(i as f64, 0))).collect()
    }

    /// Reconstruct Chern classes from Chern roots.
    /// c(E) = Π(1 + x_i) where x_i are the Chern roots.
    pub fn chern_classes_from_roots(roots: &[Polynomial]) -> Vec<Polynomial> {
        let mut total = Polynomial::one();
        for root in roots {
            let factor = Polynomial::one().add(root);
            total = total.mul(&factor);
        }

        let rank = roots.len();
        let mut classes = vec![Polynomial::one()];
        for k in 1..=rank {
            classes.push(Polynomial::monomial(total.coeff(k), k));
        }
        classes
    }

    /// Express the Todd class in terms of Chern roots.
    /// td(E) = Π x_i / (1 - e^{-x_i})
    pub fn todd_from_roots(roots: &[Polynomial]) -> Polynomial {
        let mut td = Polynomial::one();
        for root in roots {
            // Approximate: td(x) ≈ 1 + x/2 + x²/12
            let factor = Polynomial::one()
                .add(&root.scale(0.5))
                .add(&root.mul(root).scale(1.0 / 12.0));
            td = td.mul(&factor);
        }
        td
    }

    /// Verify the Whitney sum formula using splitting principle.
    pub fn verify_whitney_sum(rank_e: usize, rank_f: usize) -> bool {
        let roots_e = Self::chern_roots(rank_e);
        let roots_f = Self::chern_roots(rank_f);

        let classes_e = Self::chern_classes_from_roots(&roots_e);
        let classes_f = Self::chern_classes_from_roots(&roots_f);

        let bundle_e = VectorBundle {
            rank: rank_e,
            name: "E".to_string(),
            chern_classes: classes_e,
        };
        let bundle_f = VectorBundle {
            rank: rank_f,
            name: "F".to_string(),
            chern_classes: classes_f,
        };

        let sum = ChernClass::whitney_sum(&bundle_e, &bundle_f);
        let all_roots: Vec<_> = roots_e.into_iter().chain(roots_f).collect();
        let expected = Self::chern_classes_from_roots(&all_roots);

        // Compare total Chern classes
        let got = sum.total_chern_class();
        let mut exp_total = Polynomial::zero();
        for (_k, ck) in expected.iter().enumerate() {
            exp_total = exp_total.add(ck);
        }

        got.sub(&exp_total).coeffs.iter().all(|&c| c.abs() < 1e-10)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_zero() {
        let p = Polynomial::zero();
        assert_eq!(p.evaluate(5.0), 0.0);
        assert!(p.coeffs.is_empty());
    }

    #[test]
    fn test_polynomial_one() {
        let p = Polynomial::one();
        assert_eq!(p.evaluate(5.0), 1.0);
        assert_eq!(p.degree(), 0);
    }

    #[test]
    fn test_polynomial_add() {
        let a = Polynomial::new(vec![1.0, 2.0]); // 1 + 2t
        let b = Polynomial::new(vec![3.0, 4.0]); // 3 + 4t
        let sum = a.add(&b);
        assert_eq!(sum.coeff(0), 4.0);
        assert_eq!(sum.coeff(1), 6.0);
    }

    #[test]
    fn test_polynomial_multiply() {
        let a = Polynomial::new(vec![1.0, 1.0]); // 1 + t
        let b = Polynomial::new(vec![1.0, 1.0]); // 1 + t
        let prod = a.mul(&b);
        // (1+t)² = 1 + 2t + t²
        assert_eq!(prod.coeff(0), 1.0);
        assert_eq!(prod.coeff(1), 2.0);
        assert_eq!(prod.coeff(2), 1.0);
    }

    #[test]
    fn test_polynomial_derivative() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]); // 1 + 2t + 3t²
        let d = p.derivative();
        assert_eq!(d.coeff(0), 2.0);
        assert_eq!(d.coeff(1), 6.0);
    }

    #[test]
    fn test_trivial_bundle() {
        let b = VectorBundle::trivial(3);
        assert_eq!(b.rank, 3);
        assert_eq!(b.get_chern_class(0), Polynomial::one());
        assert_eq!(b.get_chern_class(1), Polynomial::zero());
    }

    #[test]
    fn test_line_bundle() {
        let lb = VectorBundle::line_bundle(Polynomial::monomial(5.0, 1));
        assert_eq!(lb.rank, 1);
        assert_eq!(lb.get_chern_class(1).coeff(1), 5.0);
    }

    #[test]
    fn test_whitney_sum() {
        let l1 = VectorBundle::line_bundle(Polynomial::monomial(1.0, 1));
        let l2 = VectorBundle::line_bundle(Polynomial::monomial(2.0, 1));
        let sum = ChernClass::whitney_sum(&l1, &l2);
        assert_eq!(sum.rank, 2);
        // c_1 = c_1(L1) + c_1(L2) = t + 2t = 3t
        assert_eq!(sum.get_chern_class(1).coeff(1), 3.0);
    }

    #[test]
    fn test_tangent_bundle_cp1() {
        let tb = VectorBundle::tangent_bundle_cp(1);
        assert_eq!(tb.rank, 1);
        // T(CP^1) has c_1 = 2h
        assert_eq!(tb.get_chern_class(1).coeff(1), 2.0);
    }

    #[test]
    fn test_tangent_bundle_cp2() {
        let tb = VectorBundle::tangent_bundle_cp(2);
        assert_eq!(tb.rank, 2);
        // c(T) = (1+h)^3 = 1 + 3h + 3h²
        assert_eq!(tb.get_chern_class(1).coeff(1), 3.0);
        assert_eq!(tb.get_chern_class(2).coeff(2), 3.0);
    }

    #[test]
    fn test_chern_character() {
        let lb = VectorBundle::line_bundle(Polynomial::monomial(1.0, 1));
        let ch = ChernClass::chern_character(&lb);
        // ch = 1 + c1 + c1²/2 = 1 + t + t²/2
        assert_eq!(ch.coeff(0), 1.0);
        assert_eq!(ch.coeff(1), 1.0);
        assert_eq!((ch.coeff(2) - 0.5).abs() < 1e-10, true);
    }

    #[test]
    fn test_todd_class() {
        let lb = VectorBundle::line_bundle(Polynomial::monomial(1.0, 1));
        let td = ToddClass::compute(&lb);
        // td ≈ 1 + c1/2 + c1²/12 = 1 + 0.5t + t²/12
        assert_eq!(td.coeff(0), 1.0);
        assert_eq!((td.coeff(1) - 0.5).abs() < 1e-10, true);
    }

    #[test]
    fn test_pontryagin_p1() {
        let tb = VectorBundle::tangent_bundle_cp(2);
        let p1 = PontryaginClass::p1(&tb);
        // p_1 = -c_2 = -3t²
        assert_eq!(p1.coeff(2), -3.0);
    }

    #[test]
    fn test_splitting_principle() {
        let roots = SplittingPrinciple::chern_roots(2);
        assert_eq!(roots.len(), 2);
        let classes = SplittingPrinciple::chern_classes_from_roots(&roots);
        assert_eq!(classes.len(), 3); // c_0, c_1, c_2
    }

    #[test]
    fn test_verify_whitney_sum() {
        assert!(SplittingPrinciple::verify_whitney_sum(1, 1));
    }

    #[test]
    fn test_c1_tensor_line() {
        let l1 = VectorBundle::line_bundle(Polynomial::monomial(3.0, 1));
        let l2 = VectorBundle::line_bundle(Polynomial::monomial(4.0, 1));
        let c1 = ChernClass::c1_tensor_line(&l1, &l2);
        assert_eq!(c1.coeff(1), 7.0); // 3 + 4 = 7
    }

    #[test]
    fn test_polynomial_display() {
        let p = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let s = format!("{}", p);
        assert!(s.contains("1") && s.contains("2t") && s.contains("3t^2"));
    }

    #[test]
    fn test_a_genus() {
        let tb = VectorBundle::tangent_bundle_cp(2);
        let a = PontryaginClass::a_genus(&tb);
        // Â = 1 - p_1/24 + ... = 1 + 3/24 + ...
        assert_eq!(a.coeff(0), 1.0);
        assert!((a.coeff(2) - (3.0 / 24.0)).abs() < 1e-10);
    }
}
