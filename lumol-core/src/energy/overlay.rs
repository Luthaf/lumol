// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use crate::{Potential, PairPotential, BondPotential, AnglePotential, DihedralPotential};
use crate::{Matrix3, Vector3D};

/// Add multiple [pair potentials] to create a new one as the sum of the
/// original potentials.
///
/// [pair potentials]: trait.PairPotential.html
///
/// # Examples
/// ```
/// # use lumol_core::{Potential, PairPotentialOverlay, Harmonic, LennardJones};
/// let mut overlay = PairPotentialOverlay::new();
///
/// let lj = LennardJones{ epsilon: 42.0, sigma: 2.4 };
/// overlay.add_potential(Box::new(lj.clone()));
///
/// let harmonic = Harmonic{ x0: 1.8, k: 1200.0 };
/// overlay.add_potential(Box::new(harmonic.clone()));
///
/// assert_eq!(overlay.energy(1.5), lj.energy(1.5) + harmonic.energy(1.5));
/// assert_eq!(overlay.force(1.5), lj.force(1.5) + harmonic.force(1.5));
/// ```
#[derive(Clone)]
pub struct PairPotentialOverlay {
    potentials: Vec<Box<dyn PairPotential>>,
}

impl PairPotentialOverlay {
    /// Create a new empty overlay potential.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::{Potential, PairPotentialOverlay};
    /// let overlay = PairPotentialOverlay::new();
    ///
    /// assert_eq!(overlay.energy(1.0), 0.0);
    /// ```
    pub fn new() -> PairPotentialOverlay {
        PairPotentialOverlay {
            potentials: Vec::new(),
        }
    }

    /// Add a `potential` to this overlay potential.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::{Potential, PairPotentialOverlay, Harmonic};
    /// let mut overlay = PairPotentialOverlay::new();
    ///
    /// let harmonic = Harmonic{ x0: 1.8, k: 1200.0 };
    /// overlay.add_potential(Box::new(harmonic.clone()));
    /// assert_eq!(overlay.energy(1.0), harmonic.energy(1.0));
    ///
    /// overlay.add_potential(Box::new(harmonic.clone()));
    /// assert_eq!(overlay.energy(1.0), 2.0 * harmonic.energy(1.0));
    /// ```
    pub fn add_potential(&mut self, potential: Box<dyn PairPotential>) {
        self.potentials.push(potential);
    }
}

impl Potential for PairPotentialOverlay {
    fn energy(&self, r: f64) -> f64 {
        let mut energy = 0.0;
        for potential in &self.potentials {
            energy += potential.energy(r);
        }
        return energy;
    }

    fn force(&self, r: f64) -> f64 {
        let mut force = 0.0;
        for potential in &self.potentials {
            force += potential.force(r);
        }
        return force;
    }
}

impl PairPotential for PairPotentialOverlay {
    fn virial(&self, r: &Vector3D) -> Matrix3 {
        let mut virial = Matrix3::zero();
        for potential in &self.potentials {
            virial += potential.virial(r);
        }
        return virial;
    }

    fn tail_energy(&self, cutoff: f64) -> f64 {
        let mut energy = 0.0;
        for potential in &self.potentials {
            energy += potential.tail_energy(cutoff);
        }
        return energy;
    }

    fn tail_virial(&self, cutoff: f64) -> f64 {
        let mut virial = 0.0;
        for potential in &self.potentials {
            virial += potential.tail_virial(cutoff);
        }
        return virial;
    }
}

/// Add multiple [bond potentials] to create a new one as the sum of the
/// original potentials.
///
/// [bond potentials]: trait.BondPotential.html
///
/// # Examples
/// ```
/// # use lumol_core::{Potential, BondPotentialOverlay, Harmonic, Morse};
/// let mut overlay = BondPotentialOverlay::new();
///
/// unimplemented!()
///
/// ```
#[derive(Clone)]
pub struct BondPotentialOverlay {
    potentials: Vec<Box<dyn BondPotential>>,
}

impl BondPotentialOverlay {
    /// Create a new empty overlay potential.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::{Potential, BondPotentialOverlay};
    /// let overlay = BondPotentialOverlay::new();
    ///
    /// assert_eq!(overlay.energy(1.0), 0.0);
    /// ```
    pub fn new() -> BondPotentialOverlay {
        BondPotentialOverlay {
            potentials: Vec::new(),
        }
    }

    /// Add a `potential` to this overlay potential.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::{Potential, BondPotentialOverlay, Harmonic};
    /// let mut overlay = BondPotentialOverlay::new();
    ///
    /// let harmonic = Harmonic{ x0: 1.8, k: 1200.0 };
    /// overlay.add_potential(Box::new(harmonic.clone()));
    /// assert_eq!(overlay.energy(1.0), harmonic.energy(1.0));
    ///
    /// overlay.add_potential(Box::new(harmonic.clone()));
    /// assert_eq!(overlay.energy(1.0), 2.0 * harmonic.energy(1.0));
    /// ```
    pub fn add_potential(&mut self, potential: Box<dyn BondPotential>) {
        self.potentials.push(potential);
    }
}

impl Potential for BondPotentialOverlay {
    fn energy(&self, r: f64) -> f64 {
        let mut energy = 0.0;
        for potential in &self.potentials {
            energy += potential.energy(r);
        }
        return energy;
    }

    fn force(&self, r: f64) -> f64 {
        let mut force = 0.0;
        for potential in &self.potentials {
            force += potential.force(r);
        }
        return force;
    }
}

impl BondPotential for BondPotentialOverlay {}

/// Add multiple [angle potentials] to create a new one as the sum of the
/// original potentials.
///
/// [angle potentials]: trait.AnglePotential.html
///
/// # Examples
/// ```
/// # use lumol_core::{Potential, AnglePotentialOverlay, Harmonic, Morse};
/// let mut overlay = AnglePotentialOverlay::new();
///
/// unimplemented!()
///
/// ```
#[derive(Clone)]
pub struct AnglePotentialOverlay {
    potentials: Vec<Box<dyn AnglePotential>>,
}

impl AnglePotentialOverlay {
    /// Create a new empty overlay potential.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::{Potential, AnglePotentialOverlay};
    /// let overlay = AnglePotentialOverlay::new();
    ///
    /// assert_eq!(overlay.energy(1.0), 0.0);
    /// ```
    pub fn new() -> AnglePotentialOverlay {
        AnglePotentialOverlay {
            potentials: Vec::new(),
        }
    }

    /// Add a `potential` to this overlay potential.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::{Potential, AnglePotentialOverlay, Harmonic};
    /// let mut overlay = AnglePotentialOverlay::new();
    ///
    /// let harmonic = Harmonic{ x0: 1.8, k: 1200.0 };
    /// overlay.add_potential(Box::new(harmonic.clone()));
    /// assert_eq!(overlay.energy(1.0), harmonic.energy(1.0));
    ///
    /// overlay.add_potential(Box::new(harmonic.clone()));
    /// assert_eq!(overlay.energy(1.0), 2.0 * harmonic.energy(1.0));
    /// ```
    pub fn add_potential(&mut self, potential: Box<dyn AnglePotential>) {
        self.potentials.push(potential);
    }
}

impl Potential for AnglePotentialOverlay {
    fn energy(&self, r: f64) -> f64 {
        let mut energy = 0.0;
        for potential in &self.potentials {
            energy += potential.energy(r);
        }
        return energy;
    }

    fn force(&self, r: f64) -> f64 {
        let mut force = 0.0;
        for potential in &self.potentials {
            force += potential.force(r);
        }
        return force;
    }
}

impl AnglePotential for AnglePotentialOverlay {}

/// Add multiple [dihedral potentials] to create a new one as the sum of the
/// original potentials.
///
/// [dihedral potentials]: trait.DihedralPotential.html
///
/// # Examples
/// ```
/// # use lumol_core::{Potential, DihedralPotentialOverlay, Harmonic, Morse};
/// let mut overlay = DihedralPotentialOverlay::new();
///
/// unimplemented!()
///
/// ```
#[derive(Clone)]
pub struct DihedralPotentialOverlay {
    potentials: Vec<Box<dyn DihedralPotential>>,
}

impl DihedralPotentialOverlay {
    /// Create a new empty overlay potential.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::{Potential, DihedralPotentialOverlay};
    /// let overlay = DihedralPotentialOverlay::new();
    ///
    /// assert_eq!(overlay.energy(1.0), 0.0);
    /// ```
    pub fn new() -> DihedralPotentialOverlay {
        DihedralPotentialOverlay {
            potentials: Vec::new(),
        }
    }

    /// Add a `potential` to this overlay potential.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::{Potential, AnglePotentialOverlay, Harmonic};
    /// let mut overlay = AnglePotentialOverlay::new();
    ///
    /// let harmonic = Harmonic{ x0: 1.8, k: 1200.0 };
    /// overlay.add_potential(Box::new(harmonic.clone()));
    /// assert_eq!(overlay.energy(1.0), harmonic.energy(1.0));
    ///
    /// overlay.add_potential(Box::new(harmonic.clone()));
    /// assert_eq!(overlay.energy(1.0), 2.0 * harmonic.energy(1.0));
    /// ```
    pub fn add_potential(&mut self, potential: Box<dyn DihedralPotential>) {
        self.potentials.push(potential);
    }
}

impl Potential for DihedralPotentialOverlay {
    fn energy(&self, r: f64) -> f64 {
        let mut energy = 0.0;
        for potential in &self.potentials {
            energy += potential.energy(r);
        }
        return energy;
    }

    fn force(&self, r: f64) -> f64 {
        let mut force = 0.0;
        for potential in &self.potentials {
            force += potential.force(r);
        }
        return force;
    }
}

impl DihedralPotential for DihedralPotentialOverlay {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Harmonic, LennardJones};

    #[test]
    fn pair_overlay() {
        let mut overlay = PairPotentialOverlay::new();
        assert_eq!(overlay.energy(1.0), 0.0);
        assert_eq!(overlay.force(1.0), 0.0);

        let first = Harmonic{ x0: 0.8, k: 1200.0 };
        overlay.add_potential(Box::new(first.clone()));
        assert_eq!(overlay.energy(1.0), first.energy(1.0));
        assert_eq!(overlay.force(1.0), first.force(1.0));

        let second = Harmonic{ x0: 0.9, k: 300.0 };
        overlay.add_potential(Box::new(second.clone()));
        assert_eq!(overlay.energy(1.0), first.energy(1.0) + second.energy(1.0));
        assert_eq!(overlay.force(1.0), first.force(1.0) + second.force(1.0));

        // harmonic potential has no tail correction
        assert_eq!(overlay.tail_energy(1.0), 0.0);
        assert_eq!(overlay.tail_virial(1.0), 0.0);

        let mut overlay = PairPotentialOverlay::new();
        let first = LennardJones{ sigma: 1.8, epsilon: 120.0 };
        let second = LennardJones{ sigma: 2.3, epsilon: 300.0 };
        overlay.add_potential(Box::new(first.clone()));
        overlay.add_potential(Box::new(second.clone()));

        assert_eq!(overlay.energy(1.0), first.energy(1.0) + second.energy(1.0));
        assert_eq!(overlay.force(1.0), first.force(1.0) + second.force(1.0));

        assert_eq!(overlay.tail_energy(10.0), first.tail_energy(10.0) + second.tail_energy(10.0));
        assert_eq!(overlay.tail_virial(10.0), first.tail_virial(10.0) + second.tail_virial(10.0));
    }

    #[test]
    fn bond_overlay() {
        let mut overlay = BondPotentialOverlay::new();
        assert_eq!(overlay.energy(1.0), 0.0);
        assert_eq!(overlay.force(1.0), 0.0);

        let first = Harmonic{ x0: 0.8, k: 1200.0 };
        overlay.add_potential(Box::new(first.clone()));
        assert_eq!(overlay.energy(1.0), first.energy(1.0));
        assert_eq!(overlay.force(1.0), first.force(1.0));

        let second = Harmonic{ x0: 0.9, k: 300.0 };
        overlay.add_potential(Box::new(second.clone()));
        assert_eq!(overlay.energy(1.0), first.energy(1.0) + second.energy(1.0));
        assert_eq!(overlay.force(1.0), first.force(1.0) + second.force(1.0));
    }

    #[test]
    fn angle_overlay() {
        let mut overlay = AnglePotentialOverlay::new();
        assert_eq!(overlay.energy(1.0), 0.0);
        assert_eq!(overlay.force(1.0), 0.0);

        let first = Harmonic{ x0: 0.8, k: 1200.0 };
        overlay.add_potential(Box::new(first.clone()));
        assert_eq!(overlay.energy(1.0), first.energy(1.0));
        assert_eq!(overlay.force(1.0), first.force(1.0));

        let second = Harmonic{ x0: 0.9, k: 300.0 };
        overlay.add_potential(Box::new(second.clone()));
        assert_eq!(overlay.energy(1.0), first.energy(1.0) + second.energy(1.0));
        assert_eq!(overlay.force(1.0), first.force(1.0) + second.force(1.0));
    }

    #[test]
    fn dihedral_overlay() {
        let mut overlay = DihedralPotentialOverlay::new();
        assert_eq!(overlay.energy(1.0), 0.0);
        assert_eq!(overlay.force(1.0), 0.0);

        let first = Harmonic{ x0: 0.8, k: 1200.0 };
        overlay.add_potential(Box::new(first.clone()));
        assert_eq!(overlay.energy(1.0), first.energy(1.0));
        assert_eq!(overlay.force(1.0), first.force(1.0));

        let second = Harmonic{ x0: 0.9, k: 300.0 };
        overlay.add_potential(Box::new(second.clone()));
        assert_eq!(overlay.energy(1.0), first.energy(1.0) + second.energy(1.0));
        assert_eq!(overlay.force(1.0), first.force(1.0) + second.force(1.0));
    }
}
