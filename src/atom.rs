use std::{
    io::{self, ErrorKind},
    ops::{Add, AddAssign, Neg},
    str::FromStr,
};

use approx::AbsDiffEq;
use serde::{Deserialize, Serialize};

use crate::{weights::WEIGHTS, Vec3};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Atom {
    pub atomic_number: usize,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        let eps = 1e-8;
        let close = |a: f64, b: f64| (a - b).abs() < eps;
        self.atomic_number == other.atomic_number
            && close(self.x, other.x)
            && close(self.y, other.y)
            && close(self.z, other.z)
    }
}

impl AbsDiffEq for Atom {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-8
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let close = |a: f64, b: f64| (a - b).abs() < epsilon;
        self.atomic_number == other.atomic_number
            && close(self.x, other.x)
            && close(self.y, other.y)
            && close(self.z, other.z)
    }
}

impl Neg for Atom {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            atomic_number: self.atomic_number,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Add<Vec3> for Atom {
    type Output = Atom;

    fn add(self, rhs: Vec3) -> Self::Output {
        Atom {
            x: self.x + rhs[0],
            y: self.y + rhs[1],
            z: self.z + rhs[2],
            ..self
        }
    }
}

impl AddAssign<Vec3> for Atom {
    fn add_assign(&mut self, rhs: Vec3) {
        *self = *self + rhs
    }
}

impl ToString for Atom {
    fn to_string(&self) -> String {
        format!(
            "{:2} {:15.10} {:15.10} {:15.10}",
            self.label(),
            self.x,
            self.y,
            self.z
        )
    }
}

impl FromStr for Atom {
    type Err = io::Error;

    /// parse an Atom from a line like
    ///  C 1.0 1.0 1.0
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let fields: Vec<_> = s.split_whitespace().collect();
        if fields.len() != 4 {
            return Err(io::Error::new(
                ErrorKind::Other,
                "wrong number of fields in Atom",
            ));
        }
        let coord = fields[1..].iter().map(|s| s.parse());
        if coord.clone().any(|s| s.is_err()) {
            return Err(io::Error::new(
                ErrorKind::Other,
                "failed to parse coordinate field as f64",
            ));
        }
        let coord: Vec<_> = coord.flatten().collect();
        Ok(Self::new_from_label(
            fields[0], coord[0], coord[1], coord[2],
        ))
    }
}

pub const NUMBER_TO_SYMBOL: [&str; 55] = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe",
];

impl Atom {
    pub fn new(atomic_number: usize, x: f64, y: f64, z: f64) -> Self {
        Self {
            atomic_number,
            x,
            y,
            z,
        }
    }

    pub fn new_from_label(atomic_symbol: &str, x: f64, y: f64, z: f64) -> Self {
        Self::new(
            NUMBER_TO_SYMBOL
                .iter()
                .position(|&x| x == atomic_symbol)
                .unwrap(),
            x,
            y,
            z,
        )
    }

    #[inline]
    pub const fn label(&self) -> &str {
        debug_assert!(self.atomic_number != 0 && self.atomic_number < 55);
        NUMBER_TO_SYMBOL[self.atomic_number]
    }

    pub fn coord(&self) -> Vec<f64> {
        vec![self.x, self.y, self.z]
    }

    pub fn weight(&self) -> f64 {
        WEIGHTS[self.atomic_number]
    }
}
