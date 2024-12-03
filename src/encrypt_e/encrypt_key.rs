use ark_ec::pairing::Pairing;
use std::ops::Mul;

use super::ciphertext::Ciphertext;

/// The encryption key for the encryption scheme E - ElGamal encryption.
///
/// It is a wrapper around the `bls_elgamal::EncryptKey` struct. Additionally, it
/// implements `adapt_proof` to adapt an equality proof of knowledge.
pub struct EncryptKey<E: Pairing> {
    /// The group generator.
    pub(crate) g: E::G1Affine,
    pub(crate) y1: E::G1Affine, // g^dk1
    pub(crate) y2: E::G1Affine, // g^dk2
}

impl<E: Pairing> EncryptKey<E> {
    /// Encrypt a message of (m1, m2) with randomness `v`.
    ///
    /// # Example
    ///
    /// ```rust
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use transferable_ecash::encrypt_e;
    ///
    /// type E = ark_bls12_381::Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    /// type Fr = <E as Pairing>::ScalarField;
    ///
    /// let rng = &mut test_rng();
    /// let (dk, ek) = encrypt_e::key_gen::<E, _>(rng);
    /// let (m1, m2) = (G1::rand(rng), G1::rand(rng));
    /// let v = Fr::rand(rng);
    ///
    /// let c = ek.encrypt(m1, m2, v);
    /// let (m1_, m2_) = dk.decrypt(&c);
    /// assert_eq!((m1, m2), (m1_, m2_));
    /// ```
    pub fn encrypt(&self, m1: E::G1Affine, m2: E::G1Affine, v: E::ScalarField) -> Ciphertext<E> {
        let c0 = self.g.mul(v).into();
        let c1 = (self.y1.mul(v) + m1).into();
        let c2 = (self.y2.mul(v) + m2).into();
        Ciphertext { c0, c1, c2 }
    }

    /// Randomize a ciphertext with randomness `v`.
    ///
    /// # Example
    ///
    /// ```rust
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use transferable_ecash::encrypt_e;
    ///
    /// type E = ark_bls12_381::Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    /// type Fr = <E as Pairing>::ScalarField;
    ///
    /// let rng = &mut test_rng();
    /// let (dk, ek) = encrypt_e::key_gen::<E, _>(rng);
    /// let (m1, m2) = (G1::rand(rng), G1::rand(rng));
    /// let v = Fr::rand(rng);
    /// let v_prime = Fr::rand(rng);
    ///
    /// let c = ek.encrypt(m1, m2, v);
    /// let c = ek.rerandomize(&c, v_prime);
    /// let (m1_, m2_) = dk.decrypt(&c);
    /// assert_eq!((m1, m2), (m1_, m2_));
    /// ```
    pub fn rerandomize(&self, c: &Ciphertext<E>, v: E::ScalarField) -> Ciphertext<E> {
        let c0 = (c.c0 + self.g.mul(v)).into();
        let c1 = (c.c1 + self.y1.mul(v)).into();
        let c2 = (c.c2 + self.y2.mul(v)).into();
        Ciphertext { c0, c1, c2 }
    }

    /// Verify a ciphertext with proofs.
    pub fn verify(&self) -> bool {
        todo!()
    }

    /// Adapt a proof to a rerandomization.
    ///
    /// A randomized algorithm which takes as input,
    /// - a commitment key,
    /// - an encryption public key,
    /// - a commitment,
    /// - an equality proof (i.e a Groth-Sahai proof and a commitment),
    /// - a ciphertext,
    /// - a proof,
    /// - some randomness,
    ///
    /// and outputs an equality proof.
    pub fn adapt_proof(&self) {
        todo!()
    }
}
