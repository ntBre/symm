use super::*;
use approx::assert_abs_diff_eq;
use Axis::*;

#[test]
fn test_rotate() {
    let tests = vec![
        // X around all axes
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            180.0,
            X,
        ),
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, -1.0, 0.0, 0.0)],
            180.0,
            Y,
        ),
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, -1.0, 0.0, 0.0)],
            180.0,
            Z,
        ),
        // Y around all axes
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, -1.0, 0.0)],
            180.0,
            X,
        ),
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            180.0,
            Y,
        ),
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, -1.0, 0.0)],
            180.0,
            Z,
        ),
        // Z around all axes
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, -1.0)],
            180.0,
            X,
        ),
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, -1.0)],
            180.0,
            Y,
        ),
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            180.0,
            Z,
        ),
    ];
    for test in tests {
        let h = Molecule { atoms: test.0 };
        let want = Molecule { atoms: test.1 };
        let got = h.rotate(test.2, &test.3);
        assert_eq!(got, want);
    }
}

#[test]
fn test_reflect() {
    let tests = vec![
        // X through all the planes
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, -1.0, 0.0, 0.0)],
            Plane(Y, Z),
        ),
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            Plane(X, Z),
        ),
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            Plane(X, Y),
        ),
        // Y through all the planes
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            Plane(Y, Z),
        ),
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, -1.0, 0.0)],
            Plane(X, Z),
        ),
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            Plane(X, Y),
        ),
        // Z through all the planes
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            Plane(Y, Z),
        ),
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            Plane(X, Z),
        ),
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, -1.0)],
            Plane(X, Y),
        ),
    ];
    for test in tests {
        let h = Molecule { atoms: test.0 };
        let want = Molecule { atoms: test.1 };
        let got = h.reflect(&test.2);
        assert_eq!(got, want);
    }
}

#[test]
fn test_point_group() {
    use PointGroup::*;

    struct Test<'a> {
        msg: &'a str,
        mol: &'a str,
        pg: PointGroup,
    }

    let tests = vec![
        Test {
            msg: "water",
            mol: "
  O           0.000000000    0.000000000   -0.124238453
  H           0.000000000    1.431390207    0.986041184
  H           0.000000000   -1.431390207    0.986041184
",
            pg: C2v {
                axis: Z,
                planes: vec![Plane(X, Z), Plane(Y, Z)],
            },
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
                planes: vec![Plane(Y, Z), Plane(X, Y)],
            },
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
                planes: vec![Plane(Y, Z), Plane(X, Y)],
            },
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
                axes: vec![Z, Y, X],
                planes: vec![Plane(X, Y), Plane(X, Z), Plane(Y, Z)],
            },
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
                axes: vec![Z, Y, X],
                planes: vec![Plane(X, Y), Plane(X, Z), Plane(Y, Z)],
            },
        },
    ];
    for test in tests {
        let mut mol = Molecule::from_str(test.mol).unwrap();
        mol.normalize();
        if mol.point_group() != test.pg {
            assert_eq!(
                mol.point_group(),
                test.pg,
                "wrong point group on {}",
                test.msg
            );
        }
    }
}

#[test]
fn test_irrep() {
    use Irrep::*;
    let mol_orig = Molecule::from_str(
        "
    C        0.000000   -0.888844    0.000000
    C       -0.662697    0.368254    0.000000
    C        0.662697    0.368254    0.000000
    H       -1.595193    0.906925    0.000000
    H        1.595193    0.906925    0.000000
",
    )
    .unwrap();
    let pg = mol_orig.point_group();
    let tests = vec![
        (
            vec![
                0.00, 0.03, 0.00, 0.22, -0.11, 0.00, -0.22, -0.11, 0.00, -0.57,
                0.33, 0.00, 0.57, 0.33, 0.00,
            ],
            A1,
        ),
        (
            vec![
                -0.01, 0.00, 0.00, 0.17, -0.10, 0.00, 0.17, 0.10, 0.00, -0.59,
                0.34, 0.00, -0.59, -0.34, 0.00,
            ],
            B2,
        ),
        (
            vec![
                0.00, 0.00, 0.00, 0.00, 0.00, 0.40, 0.00, 0.00, -0.40, 0.00,
                0.00, -0.58, 0.00, 0.00, 0.58,
            ],
            A2,
        ),
        (
            vec![
                0.00, 0.00, 0.16, 0.00, 0.00, -0.27, 0.00, 0.00, -0.27, 0.00,
                0.00, 0.64, 0.00, 0.00, 0.64,
            ],
            B1,
        ),
    ];
    for test in tests {
        let mol = mol_orig.clone() + test.0;
        assert_eq!(mol.irrep(&pg), test.1);
    }
}

#[test]
fn test_c2h4_irrep() {
    let mol = Molecule::from_str(
        "
C        0.0000000000       -0.0023898386        1.2600838751
H        0.0000000000        1.7483464088        2.3303799608
H        0.0000000000       -1.7425505916        2.3220592143
C        0.0000000000       -0.0014113004       -1.2600853510
H        0.0000000000        1.7444525133       -2.3255411215
H        0.0000000000       -1.7464471915       -2.3268965777
",
    )
    .unwrap();
    assert_eq!(
        mol.irrep(&PointGroup::C2v {
            axis: Z,
            planes: vec![Plane(X, Z), Plane(Y, Z)]
        }),
        Irrep::B2
    );
}

#[test]
fn test_c2h4_irrep10() {
    let mol = Molecule::from_str(
        "
C      0.00139327 -1.25400055 -0.00000000
H     -0.00036672 -2.38212447  1.70464007
H     -0.00036672 -2.38212448 -1.70464006
C     -0.00139327  1.25400055  0.00000000
H      0.00036672  2.38212448  1.70464006
H      0.00036672  2.38212448 -1.70464007
",
    )
    .unwrap();
    assert_eq!(
        mol.irrep(&PointGroup::C2v {
            axis: Y,
            planes: vec![Plane(X, Y), Plane(Y, Z)]
        }),
        Irrep::B1
    );
}

#[test]
fn test_c2h4_again() {
    let mol = Molecule::from_str(
        "C     -0.00139327 -0.00000000 -1.25400055
H      0.00036672  1.70464007 -2.38212447
H      0.00036672 -1.70464006 -2.38212448
C      0.00139327  0.00000000  1.25400055
H     -0.00036672  1.70464006  2.38212448
H     -0.00036672 -1.70464007  2.38212448
",
    )
    .unwrap();
    assert_eq!(
        mol.irrep(&PointGroup::C2v {
            axis: Z,
            planes: vec![Plane(X, Z), Plane(Y, Z)]
        }),
        Irrep::B1
    );
}

#[test]
fn test_c2h4_again_again() {
    let mol = Molecule::from_str(
        "
C     -0.00146023  0.00000002  1.36220367
H      0.00039049 -1.76535618  2.54697919
H      0.00039049  1.76535607  2.54697922
C      0.00146023 -0.00000001 -1.36220368
H     -0.00039049 -1.76535607 -2.54697922
H     -0.00039049  1.76535618 -2.54697919",
    )
    .unwrap();
    assert_eq!(
        mol.irrep_approx(
            &PointGroup::C2v {
                axis: Z,
                planes: vec![Plane(X, Z), Plane(Y, Z)]
            },
            1e-6
        )
        .unwrap(),
        Irrep::B1
    );
}
#[test]
fn test_c2h4_d2h() {
    let mol = Molecule::from_str(
        "
C      0.66679330  0.00000000 -0.42664810
H      1.23098540 -0.92361100  0.39872610
H      1.23098540  0.92361100  0.39872610
C     -0.66679330  0.00000000  0.42665160
H     -1.23098540 -0.92361100 -0.39873220
H     -1.23098540  0.92361100 -0.39873220
",
    )
    .unwrap();
    assert_eq!(
        mol.irrep_approx(
            &PointGroup::D2h {
                axes: vec![Z, X, Y],
                planes: vec![Plane(X, Y), Plane(Y, Z), Plane(X, Z)]
            },
            1e-5
        )
        .unwrap(),
        Irrep::B3g
    );
}

#[test]
fn test_com() {
    let mol = Molecule::from_str(
        "
    			H 0.0000000000 1.4313901416 0.9860410955
			O 0.0000000000 0.0000000000 -0.1242384417
			H 0.0000000000 -1.4313901416 0.9860410955
",
    )
    .unwrap();
    let got = mol.com();
    let want = Vec3::from_row_slice(&[
        0.0000000,
        0.0000000,
        9.711590454604224e-06 / 0.52917706,
    ]);
    assert_abs_diff_eq!(got, want, epsilon = 1e-8);
}

#[test]
fn test_inertia_tensor() {
    let mut mol = Molecule::from_str(
        "
    			H 0.0000000000 1.4313901416 0.9860410955
			O 0.0000000000 0.0000000000 -0.1242384417
			H 0.0000000000 -1.4313901416 0.9860410955
",
    )
    .unwrap();
    let com = mol.com();
    for atom in mol.atoms.iter_mut() {
        *atom += com;
    }
    let got = mol.inertia_tensor();
    let want = Mat3::from_row_slice(&[
        1.7743928167251328,
        0.0,
        0.0,
        0.0,
        0.61792593619882719,
        0.0,
        0.0,
        0.0,
        1.1564668805263056,
    ]) / 0.52917706
        / 0.52917706;
    assert_abs_diff_eq!(got, want, epsilon = 1e-7);
}
