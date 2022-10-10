use std::str::FromStr;

use crate::point_group::PointGroup::*;
use crate::Axis::*;
use crate::*;

struct Test<'a> {
    msg: &'a str,
    mol: &'a str,
    pg: point_group::PointGroup,
    eps: f64,
}

#[test]
fn point_group() {
    let tests = vec![
        Test {
            msg: "water",
            mol: "
  O           0.000000000    0.000000000   -0.124238453
  H           0.000000000    1.431390207    0.986041184
  H           0.000000000   -1.431390207    0.986041184
",
            pg: C2v {
                axis: Y,
                planes: [Plane(Y, Z), Plane(X, Y)],
            },
            eps: 1e-8,
        },
        Test {
            msg: "c3h2",
            mol: "
    C        0.000000   -0.888844    0.000000
    C       -0.662697    0.368254    0.000000
    C        0.662697    0.368254    0.000000
    H       -1.595193    0.906925    0.000000
    H        1.595193    0.906925    0.000000
",
            pg: C2v {
                axis: Y,
                planes: [Plane(Y, Z), Plane(X, Y)],
            },
            eps: 1e-8,
        },
        Test {
            msg: "c3h2 2",
            mol: "
         C       0.0000000000   0.0000000000   0.0000000000
         C       1.4361996439   0.0000000000   0.0000000000
         C       0.7993316223   1.1932050849   0.0000000000
         H       2.3607104536  -0.5060383602   0.0000000000
         H       0.8934572415   2.2429362063  -0.0000000000
",
            pg: C2v {
                axis: Y,
                planes: [Plane(Y, Z), Plane(X, Y)],
            },
            eps: 1e-8,
        },
        Test {
            msg: "ethylene",
            mol: "C      0.00000000 -0.00000000 -0.66360460
H     -0.00000000  0.90205573 -1.26058509
H     -0.00000000 -0.90205573 -1.26058509
C     -0.00000000  0.00000000  0.66360460
H      0.00000000  0.90205573  1.26058509
H     -0.00000000 -0.90205573  1.26058509
",
            pg: D2h {
                axes: [X, Y, Z],
                planes: [Plane(Y, Z), Plane(X, Z), Plane(X, Y)],
            },
            eps: 1e-8,
        },
        Test {
            msg: "si2c2",
            mol: "
 Si         0.0000000000        0.0000000000        1.6772191883
 Si         0.0000000000        0.0000000000       -1.6772191883
 C          0.0000000000       -0.7271526187        0.0000000000
 C          0.0000000000        0.7271526187        0.0000000000
",
            pg: D2h {
                axes: [X, Y, Z],
                planes: [Plane(Y, Z), Plane(X, Z), Plane(X, Y)],
            },
            eps: 1e-8,
        },
        Test {
            msg: "ethylene again",
            mol: "
C          0.000000000000     -0.000000000806     -0.663589004142
H          0.000000000000      0.902057584357     -1.260565246656
H         -0.000000000000     -0.902057584328     -1.260565245830
C          0.000000000000      0.000000000765      0.663589004478
H         -0.000000000000      0.902057584085      1.260565245600
H          0.000000000000     -0.902057584056      1.260565246420
",
            pg: D2h {
                axes: [X, Y, Z],
                planes: [Plane(Y, Z), Plane(X, Z), Plane(X, Y)],
            },
            eps: 1e-8,
        },
        Test {
            msg: "nh3",
            mol: "
H      0.93666628  0.31241085 -0.00000000
N      0.00001501 -0.06815209 -0.00000000
H     -0.46828832  0.31246590  0.81115089
H     -0.46828832  0.31246590 -0.81115089
",
            pg: C3v {
                axis: Z,
                // TODO two of the planes will not be aligned with Cartesian
                // planes
                plane: Plane(X, Z),
            },
            eps: 1e-6,
        },
        Test {
            msg: "c3h3+",
            mol: "
C   -0.5752253      0.5636900      0.0000000
C   -0.2005581     -0.7800038      0.0000000
C    0.7757834      0.2163138      0.0000000
H    1.8160258      0.5063683      0.0000000
H   -0.4694853     -1.8259092      0.0000000
H   -1.3465409      1.3195405      0.0000000
",
            // TODO this is actually not the right plane for C3v. I guess this
            // is why you have to rotate one of the atoms to a specific axis in
            // spectro, to make a plane coincide with the Cartesian planes. the
            // plane I'm detecting here actually indicates D3h symmetry
            pg: C3v {
                axis: Z,
                plane: Plane(X, Y),
            },
            eps: 1e-4,
        },
        // TODO this is already normalized, renormalizing messes it up
        Test {
            msg: "c3h3+ d3h",
            mol: "
        C     -0.00000736 -0.42615278  0.00000000
        C     -0.36910790  0.21308387  0.00000000
        C      0.36911526  0.21307113  0.00000000
        H      0.86400511  0.49881127  0.00000000
        H     -0.86398789  0.49884110  0.00000000
        H     -0.00001723 -0.99767867  0.00000000
        ",
            pg: D3h {
                c3: Z,
                c2: Y,
                sh: Plane(X, Y),
                sv: Plane(Y, Z),
            },
            eps: 1e-3,
        },
        // this is the pre-normalization geometry in my spectro implementation
        Test {
            msg: "bipy",
            mol: "
C      1.04734775  0.00000000  0.00000000
C      0.00000000  1.04962758  0.00000000
C      0.00000000 -0.52481379 -0.90900415
C      0.00000000 -0.52481379  0.90900415
C     -1.04734775  0.00000000  0.00000000
H      2.11983822  0.00000000  0.00000000
H     -2.11983822  0.00000000  0.00000000
",
            pg: D3h {
                c3: X,
                c2: Z,
                sh: Plane(Y, Z),
                sv: Plane(X, Z),
            },
            eps: 1e-6,
        },
    ];
    for (i, test) in tests[..].iter().enumerate() {
        let mut mol = Molecule::from_str(test.mol).unwrap();
        mol.normalize();
        let pg = mol.point_group_approx(test.eps);
        if pg != test.pg {
            println!("mol={:.8}", mol);
            assert_eq!(
                pg, test.pg,
                "wrong point group on test {i}: {}",
                test.msg
            );
        }
    }
}
