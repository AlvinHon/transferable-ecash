//! This module implements the encryption scheme E defined by  `Transferable E-cash: A Cleaner Model and the First Practical Instantiation`.

pub mod ciphertext;
pub mod decrypt_key;
pub mod encrypt_key;

use ark_ec::pairing::Pairing;
use ark_std::rand::RngCore;
use ark_std::UniformRand;

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
    let g = E::G1Affine::rand(rng);
    let dk1 = E::ScalarField::rand(rng);
    let dk2 = E::ScalarField::rand(rng);

    let dk1 = bls_elgamal::DecryptKey::new(g, dk1);
    let dk2 = bls_elgamal::DecryptKey::new(g, dk2);

    let ek1 = dk1.encrypt_key().clone();
    let ek2 = dk2.encrypt_key().clone();
    (
        DecryptKey { inner: (dk1, dk2) },
        EncryptKey { inner: (ek1, ek2) },
    )
}
