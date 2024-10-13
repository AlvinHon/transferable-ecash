use ark_ec::pairing::Pairing;

use super::ciphertext::Ciphertext;

pub struct DecryptKey<E: Pairing> {
    // decryption keys with dk1 and dk2. They must use the same generator.
    pub(crate) inner: (
        bls_elgamal::DecryptKey<E::G1>,
        bls_elgamal::DecryptKey<E::G1>,
    ),
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
    ///
    /// let rng = &mut test_rng();
    /// let (dk, ek) = encrypt_e::key_gen::<E, _>(rng);
    /// let (m1, m2) = (G1::rand(rng), G1::rand(rng));
    /// let c = ek.encrypt(rng, m1, m2);
    /// let (m1_, m2_) = dk.decrypt(&c);
    /// assert_eq!((m1, m2), (m1_, m2_));
    /// ```
    pub fn decrypt(&self, c: &Ciphertext<E>) -> (E::G1Affine, E::G1Affine) {
        let (c1, c2) = c.into();
        let m1 = self.inner.0.decrypt(c1);
        let m2 = self.inner.1.decrypt(c2);
        (m1, m2)
    }
}
