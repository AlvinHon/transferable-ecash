use ark_ec::pairing::Pairing;
use std::ops::Neg;

use super::ciphertext::Ciphertext;

pub struct DecryptKey<E: Pairing> {
    pub(crate) dk1: E::ScalarField,
    pub(crate) dk2: E::ScalarField,
}

impl<E: Pairing> DecryptKey<E> {
    /// Decrypt a ciphertext into message (m1, m2).
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
    pub fn decrypt(&self, c: &Ciphertext<E>) -> (E::G1Affine, E::G1Affine) {
        let m1 = (c.c1 + c.c0 * self.dk1.neg()).into();
        let m2 = (c.c2 + c.c0 * self.dk2.neg()).into();
        (m1, m2)
    }
}
