//! This module provides functions related to GS proof for pairing product equation.

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_std::{rand::RngCore, Zero};
use gs_ppe::{Com, Equation, Matrix, Proof, ProofSystem, Variable};

use super::CRS;

/// Commitments and proof for GS proof.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PPEProof<E: Pairing> {
    pub(crate) c: Vec<Com<E::G1>>,
    pub(crate) d: Vec<Com<E::G2>>,
    pub(crate) proof: Proof<E>,
}

impl<E: Pairing> From<ProofSystem<E>> for PPEProof<E> {
    fn from(proof: ProofSystem<E>) -> Self {
        PPEProof {
            c: proof.c,
            d: proof.d,
            proof: proof.proof,
        }
    }
}

/// Create GS proof for pairing product equation: e(A, Y) + e(X, B) = 0.
/// This function is used by encryption function in EncryptKey.
pub(crate) fn create_ppe_proof_ayxb<E: Pairing, R: RngCore>(
    rng: &mut R,
    crs: &CRS<E>,
    a: E::G1Affine,
    y: E::G2Affine,
    x: E::G1Affine,
    b: E::G2Affine,
) -> PPEProof<E> {
    // Apply:
    //      Π e(A_i, Y_i) + Π e(X_i, B_i) + ΠΠ e(X_i, Y_j)^gamma_ij = t
    // We have:
    //  n = 1, m = 1,
    //  A1 = a, B1 = b, X = x, Y = y,
    //  gamma = [[0]]

    let var_x = Variable::new(rng, x);
    let var_y = Variable::new(rng, y);
    gs_ppe::setup(
        rng,
        &crs.cks,
        &[(a, var_y)],
        &[(var_x, b)],
        &Matrix::new(&[[E::ScalarField::zero()]]),
    )
    .into()
}

/// Check GS proof for pairing product equation: e(A, Y) + e(X, B) = 0,
/// where the proof contains commitments to X and Y (i.e. have the knowledge of X and Y).
pub(crate) fn check_ppe_proof_ayxb<E: Pairing>(
    crs: &CRS<E>,
    cp: &PPEProof<E>,
    a: E::G1Affine,
    b: E::G2Affine,
) -> bool {
    Equation::new(
        vec![a],
        vec![b],
        Matrix::new(&[[E::ScalarField::zero()]]),
        PairingOutput::zero(),
    )
    .verify(&crs.cks, &cp.c, &cp.d, &cp.proof)
}

/// Create GS proof for pairing product equation: e(X1, B1) + e(X2, B2) = T.
/// This function is used by encryption function in EncryptKey.
pub(crate) fn create_ppe_proof_xbxb<E: Pairing, R: RngCore>(
    rng: &mut R,
    crs: &CRS<E>,
    x1: E::G1Affine,
    b1: E::G2Affine,
    x2: E::G1Affine,
    b2: E::G2Affine,
) -> PPEProof<E> {
    // Apply:
    //      Π e(A_i, Y_i) + Π e(X_i, B_i) + ΠΠ e(X_i, Y_j)^gamma_ij = t
    // We have:
    //  n = 0, m = 2,
    //  A1 = 0, B1 = b1, X1 = x1, X2 = x2,
    //  gamma = [[], []] i.e. dim = (2, 0)

    let var_x1 = Variable::new(rng, x1);
    let var_x2 = Variable::new(rng, x2);
    gs_ppe::setup(
        rng,
        &crs.cks,
        &[],
        &[(var_x1, b1), (var_x2, b2)],
        &Matrix::new(&[[], []]),
    )
    .into()
}

/// Check GS proof for pairing product equation: e(X1, B1) + e(X2, B2) = T.
/// where the proof contains commitments to X1 and X2 (i.e. have the knowledge of X1 and X2).
pub(crate) fn check_ppe_proof_xbxb_t<E: Pairing>(
    crs: &CRS<E>,
    cp: &PPEProof<E>,
    b1: E::G2Affine,
    b2: E::G2Affine,
    t: PairingOutput<E>,
) -> bool {
    Equation::new(vec![], vec![b1, b2], Matrix::new(&[[], []]), t)
        .verify(&crs.cks, &cp.c, &cp.d, &cp.proof)
}
