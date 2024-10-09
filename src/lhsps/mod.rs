//! This module implements the One-time linearly homomorphic structure-preserving
//! signature from Appendix B.2 of `Transferable E-cash: A Cleaner Model and the First Practical Instantiation`.

pub mod signature;
pub mod signing_key;
pub mod verifying_key;

use ark_ec::pairing::Pairing;
use ark_std::rand::RngCore;
use ark_std::UniformRand;
use std::ops::Mul;

use signing_key::SigningKey;
use verifying_key::VerifyKey;

/// Generates key pair for the one-time linearly homomorphic structure-preserving signature.
///
/// A probabilistic algorithm taking the group parameter and an integer n denoting the dimension
/// of the message to be signed. It outputs the public verification key pk and the signing key sk.
///
/// # Example
///
/// ```rust
/// use ark_std::{test_rng, UniformRand};
/// use transferable_ecash::lhsps;
///
/// let rng = &mut test_rng();
/// let (sk, pk) = lhsps::setup::<ark_bls12_381::Bls12_381, _>(rng, 5);
/// ```
pub fn setup<E: Pairing, R: RngCore>(rng: &mut R, n: usize) -> (SigningKey<E>, VerifyKey<E>) {
    let xy: Vec<(E::ScalarField, E::ScalarField)> = (0..n)
        .map(|_| (E::ScalarField::rand(rng), E::ScalarField::rand(rng)))
        .collect();

    let gz = E::G2Affine::rand(rng);
    let gr = E::G2Affine::rand(rng);

    // pk = gz^xi + gr^yi
    let pk = xy
        .iter()
        .map(|(x, y)| (gz.mul(*x) + gr.mul(*y)).into())
        .collect();

    (SigningKey { xy }, VerifyKey { gz, gr, pk })
}
