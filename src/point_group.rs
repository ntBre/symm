use std::fmt::Display;

use serde::{Deserialize, Serialize};

use crate::plane::Plane;
use crate::Axis;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PointGroup {
    C1,
    C2 {
        axis: Axis,
    },
    Cs {
        plane: Plane,
    },
    C2v {
        axis: Axis,
        planes: [Plane; 2],
    },
    C3v {
        axis: Axis,
        plane: Plane,
    },
    C5v {
        axis: Axis,
        plane: Plane,
    },
    D2h {
        axes: [Axis; 3],
        planes: [Plane; 3],
    },

    /// sh is σₕ – the plane perpendicular to the C₃ axis, while `sv` (σᵥ) is
    /// the plane that includes the C₃ axis
    D3h {
        c3: Axis,
        c2: Axis,
        sh: Plane,
        sv: Plane,
    },
}

pub enum Pg {
    C2v,
}

impl Pg {
    /// Returns `true` if the pg is [`C2v`].
    ///
    /// [`C2v`]: Pg::C2v
    #[must_use]
    pub fn is_c2v(&self) -> bool {
        matches!(self, Self::C2v)
    }
}

impl PointGroup {
    /// return the principal symmetry axis of `self`, if it has one. None
    /// otherwise
    pub fn axis(&self) -> Option<Axis> {
        match self {
            PointGroup::C1 => None,
            PointGroup::C2 { axis } => Some(*axis),
            PointGroup::Cs { plane: _ } => None,
            PointGroup::C2v { axis, planes: _ } => Some(*axis),
            PointGroup::C3v { axis, plane: _ } => Some(*axis),
            PointGroup::D2h { axes, planes: _ } => Some(axes[0]),
            PointGroup::D3h {
                c3,
                c2: _,
                sh: _,
                sv: _,
            } => Some(*c3),
            PointGroup::C5v { axis, .. } => Some(*axis),
        }
    }

    pub fn subgroup(&self, to: Pg) -> Option<Self> {
        match self {
            PointGroup::C1 => None,
            PointGroup::C2 { axis: _ } => None,
            PointGroup::Cs { plane: _ } => None,
            PointGroup::C2v { axis: _, planes: _ } => todo!(),
            PointGroup::C3v { axis: _, plane: _ } => todo!(),
            PointGroup::D2h { axes, planes } => {
                if to.is_c2v() {
                    // need the two planes containing axis
                    Some(PointGroup::C2v {
                        axis: axes[0],
                        planes: [planes[1], planes[2]],
                    })
                } else {
                    todo!()
                }
            }
            PointGroup::D3h { .. } => todo!(),
            PointGroup::C5v { .. } => todo!(),
        }
    }

    /// Returns `true` if the point group is [`C2v`].
    ///
    /// [`C2v`]: PointGroup::C2v
    #[must_use]
    pub fn is_c2v(&self) -> bool {
        matches!(self, Self::C2v { .. })
    }

    /// Returns `true` if the point group is [`D2h`].
    ///
    /// [`D2h`]: PointGroup::D2h
    #[must_use]
    pub fn is_d2h(&self) -> bool {
        matches!(self, Self::D2h { .. })
    }
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
            PointGroup::D3h { c3, c2, sh, sv } => {
                write!(f, "D3h(C3({}), C2({}), {}, {})", c3, c2, sh, sv)
            }
            PointGroup::C5v { axis, plane } => {
                write!(f, "C5v(C5({axis}), {plane})")
            }
        }
    }
}
