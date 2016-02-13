// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! This module provides some ways to initialize the velocities in a `System`
extern crate rand;
use self::rand::distributions::{Range, Normal};
use self::rand::distributions::Sample;
use self::rand::Isaac64Rng;
use self::rand::SeedableRng;

use constants::K_BOLTZMANN;
use types::Vector3D;
use simulation::{Compute, Temperature};
use super::System;

/// Scale all velocities in the `System` such that the `system` temperature
/// is `T`.
pub fn scale(system: &mut System, T: f64) {
    let instant_temperature = Temperature.compute(system);
    let factor = f64::sqrt(T / instant_temperature);
    for particle in system {
        particle.velocity = factor * particle.velocity;
    }
}

/// Random initializer for the velocities of a system.
pub trait InitVelocities {
    /// Initialize the velocities of the system.
    fn init(&mut self, system: &mut System);
    /// Set the seed of the random number generator. The default seed is 42.
    fn seed(&mut self, seed: u64);
}

/// Initialize the velocities from a Boltzman distribution.
pub struct BoltzmanVelocities {
    T: f64,
    dist: Normal,
    rng: Isaac64Rng,
}

impl BoltzmanVelocities {
    /// This `BoltzmanVelocities` initializer will initialize velocities at
    /// temperature `T`.
    pub fn new(T: f64) -> BoltzmanVelocities {
        BoltzmanVelocities{
            T: T,
            dist: Normal::new(0.0, f64::sqrt(K_BOLTZMANN * T)),
            rng: Isaac64Rng::from_seed(&[42]),
        }
    }
}

impl InitVelocities for BoltzmanVelocities {
    fn init(&mut self, system: &mut System) {
        for particle in system.iter_mut() {
            let m_inv = 1.0 / particle.mass;
            let x = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let y = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let z = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            particle.velocity = Vector3D::new(x, y, z);
        }
        scale(system, self.T);
    }

    fn seed(&mut self, seed: u64) {
        self.rng.reseed(&[seed]);
    }
}

/// Initialize the velocities from an uniform distribution.
pub struct UniformVelocities {
    T: f64,
    dist: Range<f64>,
    rng: Isaac64Rng,
}

impl UniformVelocities {
    /// This `UniformVelocities` initializer will initialize velocities at
    /// temperature `T`.
    pub fn new(T: f64) -> UniformVelocities {
        let factor = f64::sqrt(3.0*K_BOLTZMANN * T);
        UniformVelocities{
            T: T,
            dist: Range::new(-factor, factor),
            rng: Isaac64Rng::from_seed(&[42]),
        }
    }
}

impl InitVelocities for UniformVelocities {
    fn init(&mut self, system: &mut System) {
        for particle in system.iter_mut() {
            let m_inv = 1.0 / particle.mass;
            let x = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let y = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let z = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            particle.velocity = Vector3D::new(x, y, z);
        }
        scale(system, self.T);
    }

    fn seed(&mut self, seed: u64) {
        self.rng.reseed(&[seed]);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use simulation::{Compute, Temperature};

    use system::{System, Particle};

    fn testing_system() -> System {
        let mut system = System::new();
        for _ in 0..10000 {
            system.add_particle(Particle::new("F"));
        }
        return system;
    }

    #[test]
    fn init_boltzmann() {
        let mut system = testing_system();
        let mut velocities = BoltzmanVelocities::new(300.0);
        velocities.seed(1234);
        velocities.init(&mut system);
        let T = Temperature.compute(&system);
        assert_approx_eq!(T, 300.0, 1e-9);
    }

    #[test]
    fn init_uniform() {
        let mut system = testing_system();
        let mut velocities = UniformVelocities::new(300.0);
        velocities.seed(1234);
        velocities.init(&mut system);
        let T = Temperature.compute(&system);
        assert_approx_eq!(T, 300.0, 1e-9);
    }
}
