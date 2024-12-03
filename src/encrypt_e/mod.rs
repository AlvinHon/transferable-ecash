//! This module implements idea of the encryption scheme E defined by the paper
//! `Transferable E-cash: A Cleaner Model and the First Practical Instantiation`.
//!
//! The details are not obvious from the paper (Appendix B.3), but the encryption scheme is
//! seemingly a variant of ElGamal encryption. It is not clear how GS proof should be adapted
//! to this encryption scheme. The implementation here uses GS proof for providing validity
//! of the "formness" of the ciphertext, i.e. prove you know the value `r` in ciphertext `g^r`.

pub mod ciphertext;
pub mod decrypt_key;
pub mod encrypt_key;

use ark_ec::pairing::Pairing;
use ark_std::rand::RngCore;
use ark_std::UniformRand;
use std::ops::Mul;

use decrypt_key::DecryptKey;
use encrypt_key::EncryptKey;

/// Generates key pair for the encryption scheme E - ElGamal vector encryption.
///
/// # Example
///
/// ```rust
/// use ark_std::{test_rng, UniformRand};
/// use transferable_ecash::encrypt_e;
///
/// let rng = &mut test_rng();
/// let (dk, ek) = encrypt_e::key_gen::<ark_bls12_381::Bls12_381, _>(rng);
/// ```
pub fn key_gen<E: Pairing, R: RngCore>(rng: &mut R) -> (DecryptKey<E>, EncryptKey<E>) {
    let g = E::G1::rand(rng);
    let dk1 = E::ScalarField::rand(rng);
    let dk2 = E::ScalarField::rand(rng);

    let y1 = g.mul(dk1).into();
    let y2 = g.mul(dk2).into();

    (
        DecryptKey { dk1, dk2 },
        EncryptKey {
            g: g.into(),
            y1,
            y2,
        },
    )
}
