use approx::AbsDiffEq;
pub use atom::*;
pub use irrep::*;
use na::{vector, SymmetricEigen};
pub use plane::*;
pub use point_group::*;
use rotor::Rotor;
use serde::{Deserialize, Serialize};
use std::fmt::Display;

#[cfg(test)]
mod tests;

pub mod atom;
pub mod irrep;
mod mol_traits;
mod plane;
pub mod point_group;
pub mod rotor;
mod weights;

use nalgebra as na;

pub use crate::plane::Plane;

type Vec3 = na::Vector3<f64>;
type Mat3 = na::Matrix3<f64>;

static DEBUG: bool = false;

/// angbohr and weights from spectro fortran source code
pub const ANGBOHR: f64 = 0.52917706;

// TODO expand beyond cartesian axes. an alternative formulation of this is to
// align the geometry to a cartesian axis if it doesn't start like that. I think
// rotations pretty much assume you are along the cartesian axes

// restrict these to the cartesian axes for now
#[derive(
    Debug, Default, PartialEq, Eq, Copy, Clone, Serialize, Deserialize,
)]
pub enum Axis {
    X = 0,
    Y = 1,
    #[default]
    Z = 2,
}

impl Axis {
    /// return the axes not equal to `self`
    pub fn not(&self) -> (Axis, Axis) {
        use Axis::*;
        match self {
            Axis::X => (Y, Z),
            Axis::Y => (X, Z),
            Axis::Z => (X, Y),
        }
    }

    /// return the planes containing `self`
    pub fn planes(&self) -> [Plane; 2] {
        use Axis::*;
        match self {
            Axis::X => [Plane(X, Y), Plane(X, Z)],
            Axis::Y => [Plane(X, Y), Plane(Y, Z)],
            Axis::Z => [Plane(X, Z), Plane(Y, Z)],
        }
    }
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

#[macro_export]
macro_rules! molecule {
    ($($num:ident $x:literal $y:literal $z:literal$(,)?),+) => {
	$crate::Molecule::new(vec![
	    $($crate::Atom::new_from_label(stringify!($num), $x, $y, $z),)*
	    ])
    };
}

#[derive(Clone, Default, Serialize, Deserialize)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
}

fn close(a: f64, b: f64, eps: f64) -> bool {
    f64::abs(a - b) < eps
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

    /// compute the type of molecular rotor based on the moments of inertia in
    /// `moms` to the tolerance in `eps`. These tests are taken from the
    /// [Crawford Programming
    /// Projects](https://github.com/CrawfordGroup/ProgrammingProjects/blob/master/Project%2301/hints/step7-solution.md)
    pub fn rotor_type(&self, moms: &Vec3, eps: f64) -> Rotor {
        if self.atoms.len() == 2 {
            return Rotor::Diatomic;
        }
        if moms[0] < eps {
            Rotor::Linear
        } else if close(moms[0], moms[1], eps) && close(moms[1], moms[2], eps) {
            Rotor::SphericalTop
        } else if close(moms[0], moms[1], eps) && !close(moms[1], moms[2], eps)
        {
            Rotor::OblateSymmTop
        } else if !close(moms[0], moms[1], eps) && close(moms[1], moms[2], eps)
        {
            Rotor::ProlateSymmTop
        } else {
            Rotor::AsymmTop
        }
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
            .map(|a| weights::WEIGHTS[a.atomic_number])
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

    /// calls [detect_buddies] and then returns Some of the unwrapped Vec if all
    /// of the elements are initialized and None if any of them is None
    pub fn try_detect_buddies(
        &self,
        other: &Self,
        eps: f64,
    ) -> Option<Vec<usize>> {
        let tmp: Vec<_> = self
            .detect_buddies(other, eps)
            .iter()
            .cloned()
            .flatten()
            .collect();
        if tmp.len() == other.atoms.len() {
            Some(tmp)
        } else {
            None
        }
    }

    /// compute the center of mass of `self`, assuming the most abundant isotope
    /// masses
    pub fn com(&self) -> Vec3 {
        let mut sum = 0.0;
        let mut com = Vec3::zeros();
        for atom in &self.atoms {
            let w = weights::WEIGHTS[atom.atomic_number];
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
            ret[(0, 0)] += weights::WEIGHTS[*i] * (y * y + z * z);
            ret[(1, 1)] += weights::WEIGHTS[*i] * (x * x + z * z);
            ret[(2, 2)] += weights::WEIGHTS[*i] * (x * x + y * y);
            // off-diagonal
            ret[(1, 0)] -= weights::WEIGHTS[*i] * x * y;
            ret[(2, 0)] -= weights::WEIGHTS[*i] * x * z;
            ret[(2, 1)] -= weights::WEIGHTS[*i] * y * z;
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
    pub fn normalize(&mut self) -> (Vec3, Mat3, Rotor) {
        let com = self.com();
        self.translate(-com);
        let moi = self.moi();
        let (mut pr, mut axes) = symm_eigen_decomp3(moi);
        // this is what the fortran code does, rust was okay with making it INF,
        // but this plays more nicely with the math we do later
        const TOL: f64 = 1e-5;
        let rotor = self.rotor_type(&pr, TOL);
        if rotor.is_sym_top() {
            let iaxis = if close(pr[0], pr[1], TOL) {
                3
            } else if close(pr[0], pr[2], TOL) {
                2
            } else if close(pr[1], pr[2], TOL) {
                1
            } else {
                panic!("not a symmetric top: {:.8}, {}", pr, rotor);
            };

            if iaxis == 1 {
                let egr = nalgebra::matrix![
                0.0, 0.0, -1.0;
                0.0, 1.0,  0.0;
                1.0, 0.0,  0.0;
                ];
                let atemp = axes * egr;
                axes = atemp;
                pr = vector![pr[2], pr[1], pr[0]];
            } else if iaxis == 2 {
                todo!("dist.f:419");
            }
        }
        *self = self.transform(axes.transpose());
        (pr, axes, rotor)
    }

    pub fn point_group(&self) -> point_group::PointGroup {
        self.point_group_approx(1e-8)
    }

    fn is_c2(&self, ax: Axis, eps: f64) -> bool {
        let rot = self.rotate(180.0, &ax);
        rot.abs_diff_eq(self, eps)
    }

    fn is_c3(&self, ax: Axis, eps: f64) -> bool {
        let rot = self.rotate(120.0, &ax);
        rot.abs_diff_eq(self, eps)
    }

    pub fn point_group_approx(&self, eps: f64) -> point_group::PointGroup {
        use point_group::PointGroup::*;
        use Axis::*;
        let mut axes = Vec::new();
        let mut planes = Vec::new();
        for ax in [X, Y, Z] {
            // check for C2 axis
            let rot = self.rotate(180.0, &ax);
            if rot.abs_diff_eq(self, eps) {
                axes.push(ax);
            } else if DEBUG {
                eprintln!("{}", rot - self.clone());
            }
            // check for C3 axis
            let rot = self.rotate(120.0, &ax);
            if rot.abs_diff_eq(self, eps) {
                axes.push(ax);
            } else if DEBUG {
                eprintln!("{}", rot - self.clone());
            }
        }
        for plane in [Plane(X, Y), Plane(X, Z), Plane(Y, Z)] {
            let got = self.reflect(&plane);
            if got.abs_diff_eq(self, eps) {
                planes.push(plane);
            } else if DEBUG {
                eprintln!("{}", got - self.clone());
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
                    ]
                    .into_iter()
                    .filter(|p| planes.contains(p))
                    .collect()
                } else {
                    // TODO this is probably wrong too, not sure if I'm actually
                    // looking at planes in `planes` like I wasn't above
                    let ax = axes.get(0).expect(
                        "expecting at least 1 axis for more than 1 plane",
                    );
                    let mut not_it = vec![];
                    for s in axis_sum {
                        if !axes.contains(&s.0) {
                            not_it.push(s.0);
                        }
                    }
                    // want the two axes that aren't ax, sorted by their mass in
                    // axis_sum
                    vec![Plane::new(*ax, not_it[1]), Plane::new(*ax, not_it[0])]
                }
            }
            _ => panic!("impossible number of symmetry planes"),
        };
        match (axes.len(), planes.len()) {
            (0, 1) => Cs { plane: planes[0] },
            (1, 0) => C2 { axis: axes[0] },
            // NOTE should probably check the other two planes here, but I'd
            // have to rework the planes a bit. also see the note about the
            // wrong point group for c3h3+ in the tests. to handle this properly
            // I'm probably going to have to do a lot of changing
            (1, 1) => C3v {
                axis: axes[0],
                plane: planes[0],
            },
            (2, 2) => {
                let (ax, bx) = (axes[0], axes[1]);
                let (c3, c2) = if self.is_c2(ax, eps) && self.is_c3(bx, eps) {
                    (bx, ax)
                } else if self.is_c3(ax, eps) && self.is_c2(bx, eps) {
                    (ax, bx)
                } else {
                    panic!("expected one c2 and one c3 for d3h");
                };
                let (sh, sv) = if planes[0].perp() == c3 {
                    (planes[0], planes[1])
                } else if planes[1].perp() == c3 {
                    (planes[1], planes[0])
                } else {
                    panic!("expected to find a plane ⊥ c3 axis");
                };
                D3h { c3, c2, sh, sv }
            }
            (1, 2) => C2v {
                axis: axes[0],
                planes: [planes[0], planes[1]],
            },
            (3, 3) => D2h {
                axes: axes.try_into().unwrap(),
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
        let mut ret = Vec::with_capacity(self.atoms.len());
        for (i, atom) in self.to_vecs().iter().enumerate() {
            let v = mat * atom;
            ret.push(Atom::new(self.atoms[i].atomic_number, v[0], v[1], v[2]));
        }
        Self::new(ret)
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
        pg: &point_group::PointGroup,
        eps: f64,
    ) -> Result<Irrep, SymmetryError> {
        use point_group::PointGroup::*;
        use Irrep::*;
        match pg {
            C1 => Ok(A),
            C2 { axis } => {
                let new = self.rotate(180.0, axis);
                if new.abs_diff_eq(self, eps) {
                    Ok(A)
                } else if new.rotate(180.0, axis).abs_diff_eq(self, eps) {
                    Ok(B)
                } else {
                    Err(SymmetryError::new(&format!(
                        "failed to match {} for C2",
                        axis
                    )))
                }
            }
            Cs { plane } => {
                let new = self.reflect(plane);
                if new.abs_diff_eq(self, eps) {
                    Ok(Ap)
                } else if new.reflect(plane).abs_diff_eq(self, eps) {
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
                    let new = self.rotate(180.0, axis);
                    if new.abs_diff_eq(self, eps) {
                        1 // the same
                    } else if new.rotate(180.0, axis).abs_diff_eq(self, eps) {
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
            &C3v { axis: _, plane } => {
                // defer to the Cs implementation for now to satisfy summarize
                // test
                self.irrep_approx(&Cs { plane }, eps)
            }
            &D3h { c3: _, c2, sh, sv } => {
                // defer to C2v, non-abelian point groups are hard
                self.irrep_approx(
                    &C2v {
                        axis: c2,
                        planes: [sv, sh],
                    },
                    eps,
                )
            }
        }
    }

    /// calls `irrep_approx` with the value of epsilon used in `PartialEq` for
    /// `Atom` (1e-8)
    pub fn irrep(&self, pg: &point_group::PointGroup) -> Irrep {
        self.irrep_approx(pg, 1e-8).unwrap()
    }
}

/// sort the eigenvalues and eigenvectors in decreasing order by eigenvalue
pub fn eigen_sort_dec(vals: Vec3, vecs: Mat3) -> (Vec3, Mat3) {
    eigen_sort_inner(vals, vecs, true)
}

/// if `reverse` is `true`, sort the eigenvalues into decreasing order.
/// otherwise sort them in ascending order
fn eigen_sort_inner(vals: Vec3, vecs: Mat3, reverse: bool) -> (Vec3, Mat3) {
    let mut pairs: Vec<_> = vals.iter().enumerate().collect();
    if reverse {
        pairs.sort_by(|(_, a), (_, b)| b.abs().partial_cmp(&a.abs()).unwrap());
    } else {
        pairs.sort_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap());
    }
    let vec = Vec3::from_iterator(pairs.iter().map(|i| *i.1));
    let mut mat = Mat3::zeros();
    for (i, (p, _)) in pairs.iter().enumerate() {
        mat.set_column(i, &vecs.column(*p));
    }
    (vec, mat)
}

/// sort the eigenvalues and eigenvectors in ascending order by eigenvalue
pub fn eigen_sort(vals: Vec3, vecs: Mat3) -> (Vec3, Mat3) {
    eigen_sort_inner(vals, vecs, false)
}

pub fn symm_eigen_decomp3(mat: Mat3) -> (Vec3, Mat3) {
    let SymmetricEigen {
        eigenvectors: vecs,
        eigenvalues: vals,
    } = SymmetricEigen::new(mat);
    let mut pairs: Vec<_> = vals.iter().enumerate().collect();
    pairs.sort_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());
    let (_, cols) = vecs.shape();
    let mut ret = Mat3::zeros();
    for i in 0..cols {
        ret.set_column(i, &vecs.column(pairs[i].0));
    }
    (Vec3::from_iterator(pairs.iter().map(|a| *a.1)), ret)
}

#[derive(Debug, Default, PartialEq, Eq)]
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
