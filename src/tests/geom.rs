//! tests for geometrical operations like the center of mass and the moment of
//! inertia

use crate::*;
use approx::assert_abs_diff_eq;

#[test]
fn com() {
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
fn inertia_tensor() {
    let mut mol = Molecule::from_str(
        "
    			H 0.0000000000 1.4313901416 0.9860410955
			O 0.0000000000 0.0000000000 -0.1242384417
			H 0.0000000000 -1.4313901416 0.9860410955
",
    )
    .unwrap();
    let com = mol.com();
    mol.translate(-com);
    let got = mol.moi();
    let want = na::matrix![
        1.7743928167251328, 0.0, 0.0;
        0.0, 0.61792593619882719, 0.0;
        0.0, 0.0, 1.1564668805263056;
    ] / 0.52917706
        / 0.52917706;
    assert_abs_diff_eq!(got, want, epsilon = 1e-7);
}
