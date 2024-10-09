//! This module implements the Replayable-CCA encryption scheme from Appendix B.3 of
//! of `Transferable E-cash: A Cleaner Model and the First Practical Instantiation`.

use std::vec;

use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_std::ops::Mul;
use ark_std::rand::RngCore;
use ark_std::UniformRand;
use decrypt_key::DecryptKey;
use encrypt_key::EncryptKey;
use groth_sahai::{AbstractCrs, CRS};

use crate::lhsps;

pub mod ciphertext;
pub mod decrypt_key;
pub mod encrypt_key;

pub fn key_gen<E: Pairing, R: RngCore>(rng: &mut R, n: usize) -> (DecryptKey<E>, EncryptKey<E>) {
    let crs = CRS::<E>::generate_crs(rng);

    let f = E::G1Affine::rand(rng);
    let g = crs.g1_gen;

    // alpha = [alpha1, alpha2, ..., alphan]
    let alpha = (0..n)
        .map(|_| E::ScalarField::rand(rng))
        .collect::<Vec<_>>();
    // hi = g^alpha_i for i in 1..n
    let h_alpha = alpha
        .iter()
        .map(|alpha_i| g.mul(alpha_i).into())
        .collect::<Vec<_>>();
    // **
    // v1 = [f,g,1,1,...,1]
    let mut v1 = vec![f, g];
    v1.extend_from_slice(&vec![E::G1Affine::zero(); n + 1]);
    // v2 = [1,1,1,h1,h2,...,hn]
    let mut v2 = vec![E::G1Affine::zero(); 3];
    v2.extend_from_slice(&h_alpha);

    let (sk, pk) = lhsps::setup::<E, _>(rng, n + 3);
    let sigs = (sk.sign(&v1).unwrap(), sk.sign(&v2).unwrap());

    (
        DecryptKey { alpha },
        EncryptKey {
            f,
            g,
            h_alpha,
            crs,
            sigs,
        },
    )
}
