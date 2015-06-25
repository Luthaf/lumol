/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::types::Vector3D;

/// The Particle type hold basic data about a particle in the system. It is self
/// contained, so that it will be easy to send data between parrallels
/// processes.
#[derive(Clone, Debug)]
pub struct Particle {
    /// Particle name
    name: String,
    /// Particle mass
    mass: f64,
    /// Particle positions
    position: Vector3D,
    /// Particle velocity, if needed
    velocity: Vector3D,
}


impl Particle {
    pub fn new<S: Into<String>>(name: S) -> Particle {
        // TODO: get the mass here
        Particle{name: name.into(),
                 mass: 1.0,
                 position: Vector3D::new(0.0, 0.0, 0.0),
                 velocity: Vector3D::new(0.0, 0.0, 0.0)}
    }

    pub fn name<'a>(&'a self) -> &'a str {
        &self.name
    }

    pub fn set_name<S>(&mut self, name: S) where S: Into<String> {
        self.name = name.into();
    }

    pub fn mass(&self) -> f64 {
        self.mass
    }
    pub fn set_mass(&mut self, mass: f64) {
        self.mass = mass;
    }

    pub fn position<'a>(&'a self) -> &'a Vector3D {
        &self.position
    }
    pub fn set_position(&mut self, pos: Vector3D) {
        self.position = pos;
    }
    pub fn add_position(&mut self, new_pos: Vector3D) {
        self.position = self.position + new_pos;
    }

    pub fn velocity<'a>(&'a self) -> &'a Vector3D {
        &self.velocity
    }
    pub fn set_velocity(&mut self, vel: Vector3D) {
        self.velocity = vel;
    }
    pub fn add_velocity(&mut self, new_vel: Vector3D) {
        self.velocity = self.velocity + new_vel;
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ::types::*;

    #[test]
    fn set_mass() {
        let mut part = Particle::new("O");
        assert_eq!(part.name, "O");

        part.set_mass(16.0);
        assert_eq!(part.mass(), 16.0);
    }

    #[test]
    fn set_coords() {
        let mut part = Particle::new("O");
        assert_eq!(part.position(), &Vector3D::new(0.0, 0.0, 0.0));

        part.set_position(Vector3D::new(1.0, 2.0, 3.0));
        assert_eq!(part.position(), &Vector3D::new(1.0, 2.0, 3.0));

        part.set_velocity(Vector3D::new(1.0, 2.0, 3.0));
        assert_eq!(part.velocity(), &Vector3D::new(1.0, 2.0, 3.0));
    }
}
