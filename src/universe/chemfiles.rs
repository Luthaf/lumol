/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! [Chemfiles](https://github.com/chemfiles/chemfiles/) adaptators for Cymbalum.

extern crate chemfiles;
use universe::{Particle, Universe};
use universe::cells::{UnitCell, CellType};
use types::Vector3D;

use self::chemfiles::Error;

/// Convert chemfiles types to Cymbalum types
trait ToCymbalum {
    /// Output type
    type Output;
    /// Conversion function
    fn to_cymbalum(self) -> Result<Self::Output, Error>;
}

impl ToCymbalum for chemfiles::Atom {
    type Output = Particle;
    fn to_cymbalum(self) -> Result<Particle, Error> {
        let name = try!(self.name());
        let mut part = Particle::new(name);
        let mass = try!(self.mass());
        part.mass = mass as f64;
        Ok(part)
    }
}

impl ToCymbalum for chemfiles::UnitCell {
    type Output = UnitCell;
    fn to_cymbalum(self) -> Result<UnitCell, Error> {
        let cell_type = try!(self.cell_type());
        let cell = match cell_type {
            chemfiles::CellType::Infinite => UnitCell::new(),
            chemfiles::CellType::Orthorombic => {
                let (a, b, c) = try!(self.lengths());
                UnitCell::ortho(a, b, c)
            },
            chemfiles::CellType::Triclinic => {
                let (a, b, c) = try!(self.lengths());
                let (alpha, beta, gamma) = try!(self.angles());
                UnitCell::triclinic(a, b, c, alpha, beta, gamma)
            }
        };
        Ok(cell)
    }
}

impl ToCymbalum for chemfiles::Frame {
    type Output = Universe;
    fn to_cymbalum(self) -> Result<Universe, Error> {
        let cell = try!(self.cell());
        let cell = try!(cell.to_cymbalum());
        let mut universe = Universe::from_cell(cell);
        let topology = try!(self.topology());
        let positions = try!(self.positions());
        let natoms = try!(self.natoms());
        let step = try!(self.step());

        universe.set_step(step as u64);

        for i in 0..natoms {
            let atom = try!(topology.atom(i));
            let particle = try!(atom.to_cymbalum());

            universe.add_particle(particle);
            let position = Vector3D::new(
                positions[i][0] as f64,
                positions[i][1] as f64,
                positions[i][2] as f64
            );
            universe[i].position = position;
        }

        let mut bonds = try!(topology.bonds());
        while !bonds.is_empty() {
            let bond = bonds.pop().unwrap();
            if let Some(perms) = universe.add_bond(bond[0] as usize, bond[1] as usize) {
                apply_particle_permutation(&mut bonds, perms);
            }
        }
        Ok(universe)
    }
}

fn apply_particle_permutation(bonds: &mut Vec<[usize; 2]>, perms: Vec<(usize, usize)>) {
    'bonds: for bond in bonds {
        for perm in &perms {
            if bond[0] == perm.0 {
                bond[0] = perm.1;
                continue 'bonds;
            } else if bond[1] == perm.0 {
                bond[1] = perm.1;
                continue 'bonds;
            }
        }
    }
}

/// Convert a chemfiles `Frame` to an `Universe`
pub fn frame_to_universe(frame: chemfiles::Frame) -> Result<Universe, Error> {
    frame.to_cymbalum()
}

/******************************************************************************/

/// Convert Cymbalum types to chemfiles types
trait ToChemfiles {
    /// Output type
    type Output;
    /// Conversion function
    fn to_chemfiles(&self) -> Result<Self::Output, Error>;
}

impl ToChemfiles for Particle {
    type Output = chemfiles::Atom;
    fn to_chemfiles(&self) -> Result<chemfiles::Atom, Error> {
        let mut res = try!(chemfiles::Atom::new(self.name()));
        try!(res.set_mass(self.mass as f32));
        return Ok(res);
    }
}

impl ToChemfiles for UnitCell {
    type Output = chemfiles::UnitCell;
    fn to_chemfiles(&self) -> Result<chemfiles::UnitCell, Error> {
        let res = match self.celltype() {
            CellType::Infinite => {
                unimplemented!()
            }
            CellType::Orthorombic => {
                let (a, b, c) = (self.a(), self.b(), self.c());
                try!(chemfiles::UnitCell::new(a, b, c))
            },
            CellType::Triclinic => {
                let (a, b, c) = (self.a(), self.b(), self.c());
                let (alpha, beta, gamma) = (self.alpha(), self.beta(), self.gamma());
                try!(chemfiles::UnitCell::triclinic(a, b, c, alpha, beta, gamma))
            },
        };
        return Ok(res);
    }
}

impl ToChemfiles for Universe {
    type Output = chemfiles::Frame;
    fn to_chemfiles(&self) -> Result<chemfiles::Frame, Error> {

        let natoms = self.size();
        let mut frame = try!(chemfiles::Frame::new(natoms));

        try!(frame.set_step(self.step() as usize));

        let mut topology = try!(chemfiles::Topology::new());
        let mut positions = vec![[0.0f32; 3]; natoms];
        let mut velocities = vec![[0.0f32; 3]; natoms];

        for (i, p) in self.iter().enumerate() {
            let pos = p.position;
            positions[i][0] = pos.x as f32;
            positions[i][1] = pos.y as f32;
            positions[i][2] = pos.z as f32;

            let vel = p.velocity;
            velocities[i][0] = vel.x as f32;
            velocities[i][1] = vel.y as f32;
            velocities[i][2] = vel.z as f32;

            let atom = try!(p.to_chemfiles());
            try!(topology.push(&atom));
        }
        try!(frame.set_positions(positions));
        try!(frame.set_velocities(velocities));

        for molecule in self.molecules() {
            for bond in molecule.bonds() {
                try!(topology.add_bond(bond.i(), bond.j()));
            }
        }

        try!(frame.set_topology(&topology));
        // Guessing angles and dihedrals
        try!(frame.guess_topology(false));

        let cell = try!(self.cell().to_chemfiles());
        try!(frame.set_cell(&cell));
        Ok(frame)
    }
}

/// Convert an `Universe` to a chemfiles `Frame`
pub fn universe_to_frame(universe: &Universe) -> Result<chemfiles::Frame, Error> {
    universe.to_chemfiles()
}
