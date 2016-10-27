// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use energy::Potential;
use energy::{PairPotential, BondPotential, AnglePotential, DihedralPotential};

/// The `NullPotential` always returns 0 as energy and force.
///
/// It should be used to indicate that there is no potential interaction between
/// two particles kinds.
#[derive(Clone, Copy)]
pub struct NullPotential;
impl Potential for NullPotential {
    #[inline] fn energy(&self, _: f64) -> f64 {0.0}
    #[inline] fn force(&self, _: f64) -> f64 {0.0}
}

impl PairPotential for NullPotential {}
impl BondPotential for NullPotential {}
impl AnglePotential for NullPotential {}
impl DihedralPotential for NullPotential {}

/******************************************************************************/
/// Lennard-Jones potential.
///
/// The following expression is used:
/// $$ V(r) = 4 \epsilon \left[ (\sigma/r)^{12} - (\sigma/r)^{6} \right] $$
/// where $\sigma$ is a distance constant, and $\epsilon$ an energetic constant.
#[derive(Clone, Copy)]
pub struct LennardJones {
    /// Distance constant of the Lennard-Jones potential
    pub sigma: f64,
    /// Energy constant of the Lennard-Jones potential
    pub epsilon: f64,
}

impl Potential for LennardJones {
    #[inline]
    fn energy(&self, r: f64) -> f64 {
        let s6 = f64::powi(self.sigma/r, 6);
        4.0 * self.epsilon * (f64::powi(s6, 2) - s6)
    }

    #[inline]
    fn force(&self, r: f64) -> f64 {
        let s6 = f64::powi(self.sigma/r, 6);
        -24.0 * self.epsilon * (s6 - 2.0 * f64::powi(s6, 2)) / r
    }
}

impl PairPotential for LennardJones {}

/******************************************************************************/
/// Harmonic potential.
///
/// The following expression is used:
/// $$ V(x) = \frac 12 k (x - x_0)^2 $$
/// where $x_0$ is the distance equilibrium, and $k$ the elastic constant.
#[derive(Clone, Copy)]
pub struct Harmonic {
    /// Spring constant
    pub k: f64,
    /// Equilibrium value
    pub x0: f64,
}

impl Potential for Harmonic {
    #[inline]
    fn energy(&self, x: f64) -> f64 {
        let dx = x - self.x0;
        0.5 * self.k * dx * dx
    }

    #[inline]
    fn force(&self, x: f64) -> f64 {
        self.k * (self.x0 - x)
    }
}

impl PairPotential for Harmonic {}
impl BondPotential for Harmonic {}
impl AnglePotential for Harmonic {}
impl DihedralPotential for Harmonic {}

/******************************************************************************/
/// Cosine harmonic potential.
///
/// The following expression is used:
/// $$ V(r) = \frac 12 k (\cos r - \cos x_0)^2 $$
/// where $x_0$ is the distance equilibrium, and $k$ the elastic constant.
#[derive(Clone, Copy)]
pub struct CosineHarmonic {
    /// Spring constant
    k: f64,
    /// Cosine of the equilibrium value
    cos_x0: f64,
}

impl CosineHarmonic {
    /// Create a new `CosineHarmonic` potentials, with elastic constant of `k`
    /// and equilibrium value of `x0`
    pub fn new(k: f64, x0:f64) -> CosineHarmonic {
        CosineHarmonic{k: k, cos_x0: f64::cos(x0)}
    }
}

impl Potential for CosineHarmonic {
    #[inline]
    fn energy(&self, x: f64) -> f64 {
        let dr = f64::cos(x) - self.cos_x0;
        0.5 * self.k * dr * dr
    }

    #[inline]
    fn force(&self, x: f64) -> f64 {
        self.k * (f64::cos(x) - self.cos_x0) * f64::sin(x)
    }
}

impl AnglePotential for CosineHarmonic {}
impl DihedralPotential for CosineHarmonic {}

/******************************************************************************/
/// Torsion potential.
///
/// The following expression is used:
/// $$ V(r) = k(1 + \cos(n\phi_0 - \delta)) $$
/// where $k$ is the force constant, `n` the periodicity of the potential, and
/// $\delta$ the equilibrium angle.
#[derive(Clone, Copy)]
pub struct Torsion {
    /// Force constant
    pub k: f64,
    /// Equilibrium value
    pub delta: f64,
    /// Multiplicity of the potential
    pub n: usize,
}

impl Potential for Torsion {
    #[inline]
    fn energy(&self, phi: f64) -> f64 {
        let n = self.n as f64;
        let cos = f64::cos(n * phi - self.delta);
        self.k * (1.0 + cos)
    }

    #[inline]
    fn force(&self, phi: f64) -> f64 {
        let n = self.n as f64;
        let sin = f64::sin(n * phi - self.delta);
        self.k * n * sin
    }
}

impl DihedralPotential for Torsion {}

/******************************************************************************/

#[cfg(test)]
mod tests {
    use super::*;
    use energy::Potential;
    const EPS: f64 = 1e-9;

    #[test]
    fn null() {
        let null = NullPotential;
        assert_eq!(null.energy(2.0), 0.0);
        assert_eq!(null.energy(2.5), 0.0);

        assert_eq!(null.force(2.0), 0.0);
        assert_eq!(null.force(2.5), 0.0);

        let e0 = null.energy(2.0);
        let e1 = null.energy(2.0 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, null.force(2.0), EPS);
    }

    #[test]
    fn lj() {
        let lj = LennardJones{epsilon: 0.8, sigma: 2.0};
        assert_eq!(lj.energy(2.0), 0.0);
        assert_eq!(lj.energy(2.5), -0.6189584744448002);

        assert_approx_eq!(lj.force(f64::powf(2.0, 1.0/6.0) * 2.0), 0.0);
        assert_approx_eq!(lj.force(2.5), -0.95773475733504);

        let e0 = lj.energy(4.0);
        let e1 = lj.energy(4.0 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, lj.force(4.0), 1e-6);
    }

    #[test]
    fn harmonic() {
        let harm = Harmonic{k: 50.0, x0: 2.0};
        assert_eq!(harm.energy(2.0), 0.0);
        assert_eq!(harm.energy(2.5), 6.25);

        assert_eq!(harm.force(2.0), 0.0);
        assert_eq!(harm.force(2.5), -25.0);

        let e0 = harm.energy(2.1);
        let e1 = harm.energy(2.1 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, harm.force(2.1), 1e-6);
    }

    #[test]
    fn cosine_harmonic() {
        let harm = CosineHarmonic::new(50.0, 2.0);
        assert_eq!(harm.energy(2.0), 0.0);
        let dcos = f64::cos(2.5) - f64::cos(2.0);
        assert_eq!(harm.energy(2.5), 0.5 * 50.0 * dcos * dcos);

        assert_eq!(harm.force(2.0), 0.0);
        let dcos = f64::cos(2.5) - f64::cos(2.0);
        assert_eq!(harm.force(2.5), 50.0 * dcos * f64::sin(2.5));

        let e0 = harm.energy(2.3);
        let e1 = harm.energy(2.3 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, harm.force(2.3), 1e-6);
    }

    #[test]
    fn torsion() {
        let torsion = Torsion{k: 5.0, n: 3, delta: 3.0};
        assert_eq!(torsion.energy(1.0), 10.0);
        let energy = 5.0 * (1.0 + f64::cos(3.0 * 1.1 - 3.0));
        assert_eq!(torsion.energy(1.1), energy);

        assert_eq!(torsion.force(1.0), 0.0);

        let e0 = torsion.energy(4.0);
        let e1 = torsion.energy(4.0 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, torsion.force(4.0), 1e-6);
    }
}