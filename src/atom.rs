use std::ops::Neg;

#[derive(Debug, Clone, Copy)]
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

impl Atom {
    pub fn new(atomic_number: usize, x: f64, y: f64, z: f64) -> Self {
        Self {
            atomic_number,
            x,
            y,
            z,
        }
    }
}
