use std::fmt::Display;

use crate::plane::Plane;
use crate::Axis;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PointGroup {
    C1,
    C2 { axis: Axis },
    Cs { plane: Plane },
    C2v { axis: Axis, planes: [Plane; 2] },
    C3v { axis: Axis, plane: Plane },
    D2h { axes: [Axis; 3], planes: [Plane; 3] },
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
