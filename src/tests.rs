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
        mol: &'a str,
        pg: PointGroup,
    }

    let tests = vec![
        Test {
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
            mol: "
    C        0.000000   -0.888844    0.000000
    C       -0.662697    0.368254    0.000000
    C        0.662697    0.368254    0.000000
    H       -1.595193    0.906925    0.000000
    H        1.595193    0.906925    0.000000
",
            pg: C2v {
                axis: Y,
                planes: vec![Plane(X, Y), Plane(Y, Z)],
            },
        },
        Test {
            mol: "
         C       0.0000000000   0.0000000000   0.0000000000
         C       1.4361996439   0.0000000000   0.0000000000
         C       0.7993316223   1.1932050849   0.0000000000
         H       2.3607104536  -0.5060383602   0.0000000000
         H       0.8934572415   2.2429362063  -0.0000000000
",
            pg: C2v {
                axis: Y,
                planes: vec![Plane(X, Y), Plane(Y, Z)],
            },
        },
	Test {
	    mol: "C      0.00000000 -0.00000000 -0.66360460
H     -0.00000000  0.90205573 -1.26058509
H     -0.00000000 -0.90205573 -1.26058509
C     -0.00000000  0.00000000  0.66360460
H      0.00000000  0.90205573  1.26058509
H     -0.00000000 -0.90205573  1.26058509
",
	    pg: C2v {
		axis: X,
		planes: vec![Plane(X, Y), Plane(X, Z)],
	    },
	},
    ];
    for test in tests {
        let mut mol = Molecule::from_str(test.mol).unwrap();
        mol.normalize();
        assert_eq!(mol.point_group(), test.pg,);
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
            B1,
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
            B2,
        ),
    ];
    for test in tests {
        let mol = mol_orig.clone() + test.0;
        assert_eq!(mol.irrep(&pg), test.1);
    }
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
