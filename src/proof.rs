//! This module provides functions related to GS proof for internal use.

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::AffineRepr;
use ark_std::{rand::RngCore, Zero};
use groth_sahai::statement::PPE;
use groth_sahai::verifier::Verifiable;
use groth_sahai::EquType;
use groth_sahai::{
    prover::{CProof, Provable},
    Matrix, CRS,
};

/// Create GS proof for pairing product equation: e(A, Y) + e(X, B) = 0.
/// This function is used by encryption function in EncryptKey.
pub(crate) fn create_proof_ayxb<E: Pairing, R: RngCore>(
    rng: &mut R,
    crs: &CRS<E>,
    a: E::G1Affine,
    y: E::G2Affine,
    x: E::G1Affine,
    b: E::G2Affine,
) -> Result<CProof<E>, ()> {
    // Apply:
    //      Π e(A_i, Y_i) + Π e(X_i, B_i) + ΠΠ e(X_i, Y_j)^gamma_ij = t
    // We have:
    //  n = 1, m = 1,
    //  A = [a], B = [0, b], X = [0, x], Y = [y],
    //  gamma = 0
    let xvars = vec![E::G1Affine::zero(), x];
    let yvars = vec![y];
    let a_consts = vec![a];
    let b_consts = vec![E::G2Affine::zero(), b];
    let gamma: Matrix<E::ScalarField> =
        vec![vec![E::ScalarField::zero()], vec![E::ScalarField::zero()]];
    let target: PairingOutput<E> = PairingOutput::<E>::zero();
    let equ = PPE::<E> {
        a_consts,
        b_consts,
        gamma,
        target,
    };
    let proof: CProof<E> = equ.commit_and_prove(&xvars, &yvars, &crs, rng);
    equ.verify(&proof, &crs).then(|| proof).ok_or(())
}

/// Check GS proof for pairing product equation: e(A, Y) + e(X, B) = 0,
/// where the proof contains commitments to X and Y (i.e. have the knowledge of X and Y).
pub(crate) fn check_proof_ayxb<E: Pairing>(
    crs: &CRS<E>,
    cp: &CProof<E>,
    a: E::G1Affine,
    b: E::G2Affine,
) -> bool {
    if cp.equ_proofs.is_empty() {
        return false;
    }

    let equ_proof = &cp.equ_proofs[0];
    if equ_proof.equ_type != EquType::PairingProduct {
        return false;
    }

    let a_consts = vec![a];
    let b_consts = vec![E::G2Affine::zero(), b];
    let target: PairingOutput<E> = PairingOutput::<E>::zero();
    let gamma: Matrix<E::ScalarField> =
        vec![vec![E::ScalarField::zero()], vec![E::ScalarField::zero()]];

    let equ = PPE::<E> {
        a_consts,
        b_consts,
        gamma,
        target,
    };

    equ.verify(&cp, &crs)
}

/// Create GS proof for pairing product equation: e(X1, B1) + e(X2, B2) = T.
/// This function is used by encryption function in EncryptKey.
pub(crate) fn create_proof_xbxb_t<E: Pairing, R: RngCore>(
    rng: &mut R,
    crs: &CRS<E>,
    x1: E::G1Affine,
    b1: E::G2Affine,
    x2: E::G1Affine,
    b2: E::G2Affine,
    target: PairingOutput<E>,
) -> Result<CProof<E>, ()> {
    // Apply:
    //      Π e(A_i, Y_i) + Π e(X_i, B_i) + ΠΠ e(X_i, Y_j)^gamma_ij = t
    // We have:
    //  n = 0, m = 2,
    //  A = [0], B = [0, b1, b2], X = [0, x1, x2], Y = [0],
    //  gamma = 0
    let xvars = vec![E::G1Affine::zero(), x1, x2];
    let yvars = vec![E::G2Affine::zero()];
    let a_consts = vec![E::G1Affine::zero()];
    let b_consts = vec![E::G2Affine::zero(), b1, b2];
    let gamma: Matrix<E::ScalarField> = vec![
        vec![E::ScalarField::zero()],
        vec![E::ScalarField::zero()],
        vec![E::ScalarField::zero()],
    ];
    let equ = PPE::<E> {
        a_consts,
        b_consts,
        gamma,
        target,
    };
    let proof: CProof<E> = equ.commit_and_prove(&xvars, &yvars, &crs, rng);
    equ.verify(&proof, &crs).then(|| proof).ok_or(())
}

/// Check GS proof for pairing product equation: e(X1, B1) + e(X2, B2) = T.
/// where the proof contains commitments to X1 and X2 (i.e. have the knowledge of X1 and X2).
pub(crate) fn check_proof_xbxb_t<E: Pairing>(
    crs: &CRS<E>,
    cp: &CProof<E>,
    b1: E::G2Affine,
    b2: E::G2Affine,
    t: PairingOutput<E>,
) -> bool {
    if cp.equ_proofs.is_empty() {
        return false;
    }

    let equ_proof = &cp.equ_proofs[0];
    if equ_proof.equ_type != EquType::PairingProduct {
        return false;
    }

    let a_consts = vec![E::G1Affine::zero()];
    let b_consts = vec![E::G2Affine::zero(), b1, b2];
    let gamma: Matrix<E::ScalarField> = vec![
        vec![E::ScalarField::zero()],
        vec![E::ScalarField::zero()],
        vec![E::ScalarField::zero()],
    ];

    let equ = PPE::<E> {
        a_consts,
        b_consts,
        gamma,
        target: t,
    };

    equ.verify(&cp, &crs)
}
