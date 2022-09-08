use std::{
    collections::HashMap,
    fmt::Display,
    ops::{Add, BitXor},
    str::FromStr,
    string::ParseError,
};

#[cfg(test)]
mod tests;

use approx::AbsDiffEq;
pub use atom::*;
pub mod atom;
pub mod irrep;
pub use irrep::*;

use nalgebra as na;

type Vec3 = na::Vector3<f64>;
type Mat3 = na::Matrix3<f64>;

/// from <https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0>
pub const ANGBOHR: f64 = 0.5291_772_109;

/// atomic weights from https://physics.nist.gov
const WEIGHTS: [f64; 19] = [
    0.0,
    1.007_825_032,
    3.016_029_320,
    7.016_003_43,
    9.012_183_065,
    11.009_305_36,
    12.0000000,
    14.003_074_004_43,
    15.994_914_619_57,
    18.998_403_162_73,
    19.992_440_1762,
    // sodium
    22.989_769_2820,
    23.985_041_697,
    26.981_538_53,
    27.976_926_534_65,
    30.973761998,
    31.972071174,
    34.968852682,
    39.9623831237,
    // potassium
];

// TODO expand beyond cartesian axes. an alternative formulation of this is to
// align the geometry to a cartesian axis if it doesn't start like that. I think
// rotations pretty much assume you are along the cartesian axes

// restrict these to the cartesian axes for now
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum Axis {
    X = 0,
    Y = 1,
    Z = 2,
}

impl Display for Axis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Axis({})",
            match self {
                Axis::X => "X",
                Axis::Y => "Y",
                Axis::Z => "Z",
            }
        )
    }
}

impl BitXor<Axis> for Plane {
    type Output = Axis;

    fn bitxor(self, rhs: Axis) -> Self::Output {
        let Plane(ax, bx) = self;
        use Axis::*;
        match (ax, bx) {
            (X, Y) | (Y, X) => match rhs {
                X => Y,
                Y => X,
                Z => panic!("Z not in XY"),
            },
            (X, Z) | (Z, X) => match rhs {
                X => Z,
                Y => panic!("Y not in XZ"),
                Z => X,
            },
            (Y, Z) | (Z, Y) => match rhs {
                X => panic!("X not in YZ"),
                Y => Z,
                Z => Y,
            },
            _ => panic!("impossible Axis combination for Plane"),
        }
    }
}

// restrict these to combinations of cartesian axes for now. a more general
// plane is described by (a, b, c) in the equation ax + by + cz = 0
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct Plane(Axis, Axis);

impl Plane {
    /// return a normalized version of Plane
    fn new(ax: Axis, bx: Axis) -> Self {
        use Axis::*;
        match (ax, bx) {
            (X, Y) | (Y, X) => Plane(X, Y),
            (X, Z) | (Z, X) => Plane(X, Z),
            (Y, Z) | (Z, Y) => Plane(Y, Z),
            _ => panic!("impossible Axis combination for Plane"),
        }
    }

    /// return the axis perpendicular to `self`
    pub fn perp(&self) -> Axis {
        let Plane(ax, bx) = self;
        use Axis::*;
        match (ax, bx) {
            (X, Y) | (Y, X) => Z,
            (X, Z) | (Z, X) => Y,
            (Y, Z) | (Z, Y) => X,
            _ => panic!("impossible Axis combination for Plane"),
        }
    }
}

impl Display for Plane {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Plane({}, {})", self.0, self.1)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum PointGroup {
    C1,
    C2 { axis: Axis },
    Cs { plane: Plane },
    C2v { axis: Axis, planes: [Plane; 2] },
    C3v { axis: Axis, plane: Plane },
    D2h { axes: Vec<Axis>, planes: [Plane; 3] },
}

impl Display for PointGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PointGroup::C1 => write!(f, "C1"),
            PointGroup::C2 { axis: a } => write!(f, "C2({})", a),
            PointGroup::Cs { plane: p } => write!(f, "Cs({})", p),
            PointGroup::C2v {
                axis: a,
                planes: ps,
            } => write!(f, "C2v({}, {}, {})", a, ps[0], ps[1]),
            PointGroup::C3v { axis: a, plane } => {
                write!(f, "C3v({}, {})", a, plane)
            }
            PointGroup::D2h { axes, planes } => {
                write!(
                    f,
                    "D2h({}, {}, {}, {}, {}, {})",
                    axes[0], axes[1], axes[2], planes[0], planes[1], planes[2]
                )
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
}

/// A Molecule is AbsDiffEq if each of its Atoms is
impl AbsDiffEq for Molecule {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-8
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        assert!(self.atoms.len() == other.atoms.len());
        let mut theirs = other.atoms.clone();
        if self.atoms.len() != theirs.len() {
            return false;
        }
        for atom in &self.atoms {
            let mut pops = Vec::new();
            let mut found = false;
            for (i, btom) in theirs.iter().enumerate() {
                if atom.abs_diff_eq(btom, epsilon) {
                    pops.push(i);
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
            // remove high indices first
            pops.sort();
            pops.reverse();
            for p in pops {
                theirs.remove(p);
            }
        }
        true
    }
}

impl PartialEq for Molecule {
    /// compare molecules for equality, irrespective of order. try to find an
    /// atom in other that equals the current atom in self. If found, remove it,
    /// so it can't be double-counted.
    fn eq(&self, other: &Self) -> bool {
        let mut theirs = other.atoms.clone();
        if self.atoms.len() != theirs.len() {
            return false;
        }
        for atom in &self.atoms {
            let mut pops = Vec::new();
            let mut found = false;
            for (i, btom) in theirs.iter().enumerate() {
                if *atom == *btom {
                    pops.push(i);
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
            // remove high indices first
            pops.sort();
            pops.reverse();
            for p in pops {
                theirs.remove(p);
            }
        }
        true
    }
}

impl FromStr for Molecule {
    type Err = ParseError;

    /// parse lines like
    ///      O           0.000000000    0.000000000   -0.124238453
    ///      H           0.000000000    1.431390207    0.986041184
    ///      H           0.000000000   -1.431390207    0.986041184
    /// into a molecule
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut ret = Self::default();
        let atomic_symbols: HashMap<_, _> = NUMBER_TO_SYMBOL
            .iter()
            .enumerate()
            .map(|(i, s)| (s.to_string(), i))
            .collect();
        for line in s.lines() {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() == 4 {
                let sym = if let Some(&s) = atomic_symbols.get(fields[0]) {
                    s
                } else {
                    panic!(
                        "atomic symbol '{}' not found, tell Brent!",
                        fields[0]
                    );
                };
                ret.atoms.push(Atom::new(
                    sym,
                    fields[1].parse().unwrap(),
                    fields[2].parse().unwrap(),
                    fields[3].parse().unwrap(),
                ));
            }
        }
        Ok(ret)
    }
}

impl Add<Vec<f64>> for Molecule {
    type Output = Self;

    /// panics if the size of `rhs` doesn't align with the size of `self.atoms`
    fn add(self, rhs: Vec<f64>) -> Self::Output {
        if 3 * self.atoms.len() != rhs.len() {
            panic!(
                "{} atoms but {} displacements",
                self.atoms.len(),
                rhs.len()
            );
        }
        let mut ret = self.clone();
        // panic above ensures rhs is exactly divisble by 3
        for (i, chunk) in rhs.chunks_exact(3).enumerate() {
            ret.atoms[i].x += chunk[0];
            ret.atoms[i].y += chunk[1];
            ret.atoms[i].z += chunk[2];
        }
        ret
    }
}

impl Display for Molecule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let precision = f.precision().unwrap_or(8);
        let width = f.width().unwrap_or(precision + 4);
        writeln!(f)?;
        for atom in &self.atoms {
            writeln!(
                f,
                "{:5}{:w$.p$}{:w$.p$}{:w$.p$}",
                NUMBER_TO_SYMBOL[atom.atomic_number],
                atom.x,
                atom.y,
                atom.z,
                w = width,
                p = precision,
            )?;
        }
        Ok(())
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self {
            atoms: Default::default(),
        }
    }
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>) -> Self {
        Self { atoms }
    }

    fn to_vecs(&self) -> Vec<Vec3> {
        let mut ret = Vec::with_capacity(self.atoms.len());
        for atom in &self.atoms {
            ret.push(Vec3::new(atom.x, atom.y, atom.z));
        }
        ret
    }

    /// build a `Molecule` from a slice of coordinates and a slice of
    /// atomic_numbers
    pub fn from_slices(atomic_numbers: Vec<usize>, coords: &[f64]) -> Self {
        assert_eq!(3 * atomic_numbers.len(), coords.len());
        let mut atoms = Vec::new();
        for (i, atom) in coords.chunks(3).enumerate() {
            atoms.push(Atom::new(atomic_numbers[i], atom[0], atom[1], atom[2]));
        }
        Self { atoms }
    }

    /// return the atomic numbers of each atoms as a vector
    pub fn atomic_numbers(&self) -> Vec<usize> {
        self.atoms.iter().map(|a| a.atomic_number).collect()
    }

    /// return the atomic numbers of each atoms as a vector
    pub fn weights(&self) -> Vec<f64> {
        self.atoms
            .iter()
            .map(|a| WEIGHTS[a.atomic_number])
            .collect()
    }

    /// convert the coordinates in `self` from Angstroms to Bohr
    pub fn to_bohr(&mut self) {
        for atom in self.atoms.iter_mut() {
            atom.x /= ANGBOHR;
            atom.y /= ANGBOHR;
            atom.z /= ANGBOHR;
        }
    }

    /// convert the coordinates in `self` from Bohr to Angstroms
    pub fn to_angstrom(&mut self) {
        for atom in self.atoms.iter_mut() {
            atom.x *= ANGBOHR;
            atom.y *= ANGBOHR;
            atom.z *= ANGBOHR;
        }
    }

    /// a buddy is a mapping from one Vec<Atom> to another. The returned Vec
    /// contains the indices in `self` that correspond to atoms in `other`.
    /// `eps` is used in the `Atom` `AbsDiffEq` call to check the equality of
    /// two atoms.
    pub fn detect_buddies(&self, other: &Self, eps: f64) -> Vec<Option<usize>> {
        assert_eq!(self.atoms.len(), other.atoms.len());
        let mut ret = vec![None; self.atoms.len()];
        for (i, atom) in self.atoms.iter().enumerate() {
            for (j, btom) in other.atoms.iter().enumerate() {
                if atom.abs_diff_eq(btom, eps) {
                    ret[i] = Some(j);
                    break;
                }
            }
        }
        ret
    }

    /// compute the center of mass of `self`, assuming the most abundant isotope
    /// masses
    pub fn com(&self) -> Vec3 {
        let mut sum = 0.0;
        let mut com = Vec3::zeros();
        for atom in &self.atoms {
            let w = WEIGHTS[atom.atomic_number];
            sum += w;
            com += w * Vec3::from_row_slice(&atom.coord());
        }
        com / sum
    }

    /// compute the moment of inertia tensor
    pub fn moi(&self) -> Mat3 {
        let mut ret = Mat3::zeros();
        for atom in &self.atoms {
            let Atom {
                x,
                y,
                z,
                atomic_number: i,
            } = atom;
            // diagonal
            ret[(0, 0)] += WEIGHTS[*i] * (y * y + z * z);
            ret[(1, 1)] += WEIGHTS[*i] * (x * x + z * z);
            ret[(2, 2)] += WEIGHTS[*i] * (x * x + y * y);
            // off-diagonal
            ret[(1, 0)] -= WEIGHTS[*i] * x * y;
            ret[(2, 0)] -= WEIGHTS[*i] * x * z;
            ret[(2, 1)] -= WEIGHTS[*i] * y * z;
        }
        ret
    }

    /// eigenfactorize the moment of inertia tensor and return the principal
    /// moments of inertia as a Vec3
    pub fn principal_moments(&self) -> Vec3 {
        let it = self.moi();
        let sym = na::SymmetricEigen::new(it);
        sym.eigenvalues
    }

    /// eigenfactorize the moment of inertia tensor and return the principal
    /// axes as a 3x3 matrix
    pub fn principal_axes(&self) -> Mat3 {
        let it = self.moi();
        let sym = na::SymmetricEigen::new(it);
        sym.eigenvectors
    }

    /// translate each of the atoms in `self` by vec
    pub fn translate(&mut self, vec: Vec3) -> &mut Self {
        for atom in self.atoms.iter_mut() {
            *atom += vec;
        }
        self
    }

    /// normalize `self` by translating to the center of mass and orienting the
    /// molecule such that the rotational axes are aligned with the Cartesian
    /// axes. additionally, sort the axes such that x corresponds to the
    /// smallest moment of inertia and z to the largest. returns the vector of
    /// principal moments of inertia and the principal axes used to make the
    /// transformation. adapted from SPECTRO
    pub fn normalize(&mut self) -> (Vec3, Mat3) {
        let com = self.com();
        self.translate(-com);
        let pr = self.principal_moments();
        let axes = self.principal_axes();
        let (pr, axes) = eigen_sort(pr, axes);
        *self = self.transform(axes.transpose());
        (pr, axes)
    }

    pub fn point_group(&self) -> PointGroup {
        self.point_group_approx(1e-8)
    }

    pub fn point_group_approx(&self, eps: f64) -> PointGroup {
        use Axis::*;
        use PointGroup::*;
        let mut axes = Vec::new();
        let mut planes = Vec::new();
        for ax in vec![X, Y, Z] {
            // check for C2 axis
            if self.rotate(180.0, &ax).abs_diff_eq(self, eps) {
                axes.push(ax);
            }
            // check for C3 axis
            let got = self.rotate(120.0, &ax);
            if got.abs_diff_eq(self, eps) {
                axes.push(ax);
            }
        }
        for plane in vec![Plane(X, Y), Plane(X, Z), Plane(Y, Z)] {
            if self.reflect(&plane).abs_diff_eq(self, eps) {
                planes.push(plane);
            }
        }

        let axis_sum = self.axis_sum();
        let axes = if axes.len() > 1 {
            // multiple choices, so take highest mass-weighted coordinate as the
            // principal axis
            let mut axes_new = vec![];
            for s in axis_sum {
                if axes.contains(&s.0) {
                    axes_new.push(s.0);
                }
            }
            axes_new
        } else {
            axes
        };

        // if there is only one axis, you can't blindly take the highest mass
        // axis because it might not be an actual symmetry element... this only
        // works for high-symmetry groups like D2h with multiple axes AND
        // multiple planes. it doesn't work for a simple C2v molecule like
        // allyl+. this might still be wrong for even more specific cases, but
        // it works better than before
        let planes = match planes.len() {
            0 | 1 => planes,
            2 | 3 => {
                if axes.len() > 1 {
                    // the order of planes is
                    // 1. Plane(highest axis, third-highest axis)
                    // 2. Plane(highest axis, second-highest axis)
                    // 3. Plane(second-highest, third-highest)
                    [
                        Plane::new(axis_sum[0].0, axis_sum[2].0),
                        Plane::new(axis_sum[0].0, axis_sum[1].0),
                        Plane::new(axis_sum[1].0, axis_sum[2].0),
                    ][..planes.len()]
                        .to_vec()
                } else {
                    let ax = axes.get(0).expect(
                        "expecting at least 1 axis for more than 1 plane",
                    );
                    let mut not_it = vec![];
                    for s in axis_sum {
                        if !axes.contains(&s.0) {
                            not_it.push(s.0);
                        }
                    }
                    // want the two axes that aren't ax, sorted by their mass in axis_sum
                    vec![Plane::new(*ax, not_it[1]), Plane::new(*ax, not_it[0])]
                }
            }
            _ => panic!("impossible number of symmetry planes"),
        };
        match (axes.len(), planes.len()) {
            (0, 1) => Cs { plane: planes[0] },
            (1, 0) => C2 { axis: axes[0] },
            // NOTE should probably check the other two planes here, but I'd
            // have to rework the planes a bit
            (1, 1) => C3v {
                axis: axes[0],
                plane: planes[0],
            },
            (1, 2) => C2v {
                axis: axes[0],
                planes: [planes[0], planes[1]],
            },
            (3, 3) => D2h {
                axes,
                // for some reason you put the least mass plane first
                planes: [planes[2], planes[0], planes[1]],
            },
            _ => C1,
        }
    }

    /// compute the mass-weighted sum of the axes in `self` and sort them such
    /// that axis with the highest sum is first. This axis is the principal
    /// axis; combining the principal axis with the second-highest sum gives the
    /// main plane
    fn axis_sum(&self) -> [(Axis, f64); 3] {
        use Axis::*;
        let mut sum = [(X, 0.0), (Y, 0.0), (Z, 0.0)];
        for atom in &self.atoms {
            sum[0].1 += atom.weight() * atom.x.abs();
            sum[1].1 += atom.weight() * atom.y.abs();
            sum[2].1 += atom.weight() * atom.z.abs();
        }
        // reverse sort by float part
        sum.sort_by(|f, g| g.1.partial_cmp(&f.1).unwrap());
        sum
    }

    /// apply the transformation matrix `mat` to the atoms in `self` and return
    /// the new Molecule
    pub fn transform(&self, mat: na::Matrix3<f64>) -> Self {
        let mut ret = Self::default();
        for (i, atom) in self.to_vecs().iter().enumerate() {
            let v = mat * atom;
            ret.atoms.push(Atom::new(
                self.atoms[i].atomic_number,
                v[0],
                v[1],
                v[2],
            ));
        }
        ret
    }

    pub fn rotate(&self, deg: f64, axis: &Axis) -> Self {
        use Axis::*;
        let deg = deg.to_radians();
        let ct = deg.cos();
        let st = deg.sin();
        // from
        // https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
        #[rustfmt::skip]
	let rot_mat = match axis {
            X => {
		na::Matrix3::new(
		    1., 0., 0.,
		    0., ct, -st,
		    0., st, ct,
		)
            }
            Y => {
		na::Matrix3::new(
		    ct, 0., st,
		    0., 1., 0.,
		    -st, 0., ct,
		)
            }
            Z => {
		na::Matrix3::new(
		    ct, -st, 0.,
		    st, ct, 0.,
		    0., 0., 1.,
		)
            }
        };
        self.transform(rot_mat)
    }

    #[rustfmt::skip]
    /// return the special case of the Householder reflection in 3 dimensions
    /// described here:
    /// <https://en.wikipedia.org/wiki/Transformation_matrix#Reflection_2>
    fn householder(a: f64, b: f64, c: f64) -> na::Matrix3<f64> {
        na::Matrix3::new(
            1. - 2. * a * a, -2. * a * b, -2. * a * c,
            -2. * a * b, 1. - 2. * b * b, -2. * b * c,
            -2. * a * c, -2. * b * c, 1. - 2. * c * c,
        )
    }

    pub fn reflect(&self, plane: &Plane) -> Self {
        use Axis::*;
        let ref_mat = match plane {
            Plane(X, Y) | Plane(Y, X) => Self::householder(0.0, 0.0, 1.0),
            Plane(X, Z) | Plane(Z, X) => Self::householder(0.0, 1.0, 0.0),
            Plane(Y, Z) | Plane(Z, Y) => Self::householder(1.0, 0.0, 0.0),
            _ => panic!("unrecognized plane {:?}", plane),
        };
        self.transform(ref_mat)
    }

    pub fn irrep_approx(
        &self,
        pg: &PointGroup,
        eps: f64,
    ) -> Result<Irrep, SymmetryError> {
        use Irrep::*;
        use PointGroup::*;
        match pg {
            C1 => Ok(A),
            C2 { axis } => {
                let new = self.rotate(180.0, &axis);
                if new.abs_diff_eq(self, eps) {
                    Ok(A)
                } else if new.rotate(180.0, &axis).abs_diff_eq(self, eps) {
                    Ok(B)
                } else {
                    Err(SymmetryError::new(&format!(
                        "failed to match {} for C2",
                        axis
                    )))
                }
            }
            Cs { plane } => {
                let new = self.reflect(&plane);
                if new.abs_diff_eq(self, eps) {
                    Ok(Ap)
                } else if new.reflect(&plane).abs_diff_eq(self, eps) {
                    Ok(App)
                } else {
                    Err(SymmetryError::new(&format!(
                        "failed to match {plane} for Cs"
                    )))
                }
            }
            // TODO this is where the plane order can matter - B1 vs B2. as long
            // as you call `irrep` multiple times with the same PointGroup, you
            // should get consistent results at least. source of the issue is in
            // point_group - order of planes there should be based on something
            // besides random choice of iteration order in the implementation
            // (mass?)
            C2v { axis, planes } => {
                let mut chars = (0, 0, 0);
                // TODO would be nice to abstract these into some kind of apply
                // function
                chars.0 = {
                    let new = self.rotate(180.0, &axis);
                    if new.abs_diff_eq(self, eps) {
                        1 // the same
                    } else if new.rotate(180.0, &axis).abs_diff_eq(self, eps) {
                        -1 // the opposite
                    } else {
                        0 // something else
                    }
                };
                chars.1 = {
                    let new = self.reflect(&planes[0]);
                    if new.abs_diff_eq(self, eps) {
                        1
                    } else if new.reflect(&planes[0]).abs_diff_eq(self, eps) {
                        -1
                    } else {
                        0
                    }
                };
                chars.2 = {
                    let new = self.reflect(&planes[1]);
                    if new.abs_diff_eq(self, eps) {
                        1
                    } else if new.reflect(&planes[1]).abs_diff_eq(self, eps) {
                        -1
                    } else {
                        0
                    }
                };
                match chars {
                    (1, 1, 1) => Ok(A1),
                    (1, -1, -1) => Ok(A2),
                    (-1, 1, -1) => Ok(B1),
                    (-1, -1, 1) => Ok(B2),
                    _ => Err(SymmetryError::new(&format!(
                        "failed to match {:?}",
                        chars
                    ))),
                }
            }
            D2h { axes, planes } => {
                // NOTE: skipping the inversion for now

                // C2, C2, C2, σ, σ, σ
                let mut chars = (0, 0, 0, 0, 0, 0);
                // first axis
                chars.0 = {
                    let new = self.rotate(180.0, &axes[0]);
                    if new.abs_diff_eq(self, eps) {
                        1 // the same
                    } else if new.rotate(180.0, &axes[0]).abs_diff_eq(self, eps)
                    {
                        -1 // the opposite
                    } else {
                        0 // something else
                    }
                };
                // second axis
                chars.1 = {
                    let new = self.rotate(180.0, &axes[1]);
                    if new.abs_diff_eq(self, eps) {
                        1 // the same
                    } else if new.rotate(180.0, &axes[1]).abs_diff_eq(self, eps)
                    {
                        -1 // the opposite
                    } else {
                        0 // something else
                    }
                };
                // third axis
                chars.2 = {
                    let new = self.rotate(180.0, &axes[2]);
                    if new.abs_diff_eq(self, eps) {
                        1 // the same
                    } else if new.rotate(180.0, &axes[2]).abs_diff_eq(self, eps)
                    {
                        -1 // the opposite
                    } else {
                        0 // something else
                    }
                };
                // first plane
                chars.3 = {
                    let new = self.reflect(&planes[0]);
                    if new.abs_diff_eq(self, eps) {
                        1
                    } else if new.reflect(&planes[0]).abs_diff_eq(self, eps) {
                        -1
                    } else {
                        0
                    }
                };
                // second plane
                chars.4 = {
                    let new = self.reflect(&planes[1]);
                    if new.abs_diff_eq(self, eps) {
                        1
                    } else if new.reflect(&planes[1]).abs_diff_eq(self, eps) {
                        -1
                    } else {
                        0
                    }
                };
                // third plane
                chars.5 = {
                    let new = self.reflect(&planes[2]);
                    if new.abs_diff_eq(self, eps) {
                        1
                    } else if new.reflect(&planes[2]).abs_diff_eq(self, eps) {
                        -1
                    } else {
                        0
                    }
                };
                match chars {
                    (1, 1, 1, 1, 1, 1) => Ok(Ag),
                    (1, -1, -1, 1, -1, -1) => Ok(B1g),
                    (-1, 1, -1, -1, 1, -1) => Ok(B2g),
                    (-1, -1, 1, -1, -1, 1) => Ok(B3g),
                    (1, 1, 1, -1, -1, -1) => Ok(Au),
                    (1, -1, -1, -1, 1, 1) => Ok(B1u),
                    (-1, 1, -1, 1, -1, 1) => Ok(B2u),
                    (-1, -1, 1, 1, 1, -1) => Ok(B3u),
                    _ => Err(SymmetryError::new(&format!(
                        "failed to match {:?} on\n{}",
                        chars, &self
                    ))),
                }
            }
            C3v { axis: _, plane: _ } => todo!(),
        }
    }

    /// calls `irrep_approx` with the value of epsilon used in `PartialEq` for
    /// `Atom` (1e-8)
    pub fn irrep(&self, pg: &PointGroup) -> Irrep {
        self.irrep_approx(pg, 1e-8).unwrap()
    }
}

/// sort the eigenvalues and eigenvectors in decreasing order by eigenvalue
pub fn eigen_sort_dec(vals: Vec3, vecs: Mat3) -> (Vec3, Mat3) {
    let mut pairs: Vec<_> = vals.iter().enumerate().collect();
    pairs.sort_by(|(_, a), (_, b)| b.partial_cmp(&a).unwrap());
    let vec = Vec3::from_iterator(pairs.iter().map(|i| i.1.clone()));
    let mut mat = Mat3::zeros();
    for i in 0..3 {
        mat.set_column(i, &vecs.column(pairs[i].0));
    }
    (vec, mat)
}

pub fn eigen_sort(vals: Vec3, vecs: Mat3) -> (Vec3, Mat3) {
    let mut pairs: Vec<_> = vals.iter().enumerate().collect();
    pairs.sort_by(|(_, a), (_, b)| a.partial_cmp(&b).unwrap());
    let vec = Vec3::from_iterator(pairs.iter().map(|i| i.1.clone()));
    let mut mat = Mat3::zeros();
    for i in 0..3 {
        mat.set_column(i, &vecs.column(pairs[i].0));
    }
    (vec, mat)
}

#[derive(Debug, Default, PartialEq)]
pub struct SymmetryError {
    msg: String,
}

impl SymmetryError {
    pub fn new(msg: &str) -> Self {
        Self {
            msg: String::from(msg),
        }
    }

    pub fn msg(&self) -> &str {
        self.msg.as_ref()
    }
}
