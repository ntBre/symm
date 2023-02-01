use std::str::FromStr;

use crate::point_group::PointGroup::*;
use crate::Axis::*;
use crate::*;

struct Test<'a> {
    msg: &'a str,
    mol: &'a str,
    pg: point_group::PointGroup,
    eps: f64,
    norm: bool,
}

#[test]
fn point_group() {
    let tests = vec![
        Test {
            msg: "cp-",
            mol: "
C     -0.07340837  1.20099830 -0.00000000
C      1.11953278  0.44094426 -0.00000000
C      0.76531752 -0.92847953  0.00000000
C     -0.64654028 -1.01477625  0.00000000
C     -1.16490165  0.30131323 -0.00000000
H     -2.19491771  0.56773722  0.00000000
H     -1.21821708 -1.91205038 -0.00000000
H      1.44201819 -1.74944935  0.00000000
H      2.10943330  0.83083121 -0.00000000
H     -0.13831669  2.26293122  0.00000000
",
            pg: C5v {
                axis: Z,
                plane: Plane(X, Y),
            },
            eps: 1e-4,
            norm: true,
        },
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
            norm: true,
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
            norm: true,
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
            norm: true,
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
                planes: [Plane(X, Y), Plane(X, Z), Plane(Y, Z)],
            },
            eps: 1e-8,
            norm: true,
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
                planes: [Plane(X, Y), Plane(X, Z), Plane(Y, Z)],
            },
            eps: 1e-8,
            norm: true,
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
                planes: [Plane(X, Y), Plane(X, Z), Plane(Y, Z)],
            },
            eps: 1e-8,
            norm: true,
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
            norm: true,
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
            norm: true,
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
            norm: true,
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
                c3: Z,
                c2: X,
                sh: Plane(X, Y),
                sv: Plane(X, Z),
            },
            eps: 1e-6,
            norm: true,
        },
        Test {
            msg: "naphthalene",
            mol: "
C     -1.24985877 -1.40744098 -0.00004085
C     -2.43325152 -0.71298818 -0.00008992
C     -2.43324132  0.71299935  0.00001093
C     -1.24984194  1.40744370  0.00008470
C      0.00000616  0.70960162  0.00004712
C      1.24986006  1.40744164 -0.00004505
C      2.43325079  0.71298767 -0.00008695
C      2.43323985 -0.71300059  0.00000945
C      1.24984291 -1.40744315  0.00008107
C     -0.00000704 -0.70960102  0.00004510
H      1.23985505 -2.49618789  0.00019554
H      3.38837675 -1.23551703  0.00004241
H      3.38840493  1.23547349 -0.00021344
H      1.23987806  2.49618830 -0.00010622
H     -1.23984739  2.49619016  0.00017846
H     -3.38837934  1.23551696  0.00004670
H     -3.38840615 -1.23547642 -0.00023892
H     -1.23987206 -2.49618816 -0.00009043
",
            pg: D2h {
                axes: [X, Y, Z],
                planes: [Plane(X, Y), Plane(X, Z), Plane(Y, Z)],
            },
            eps: 1e-3,
            norm: false,
        },
    ];
    for (i, test) in tests[..].iter().enumerate() {
        let mut mol = Molecule::from_str(test.mol).unwrap();
        if test.norm {
            mol.normalize();
        }
        let pg = mol.point_group_approx(test.eps);
        if pg != test.pg {
            println!("mol={mol:.8}");
            assert_eq!(
                pg, test.pg,
                "wrong point group on test {i}: {}",
                test.msg
            );
        }
    }
}
