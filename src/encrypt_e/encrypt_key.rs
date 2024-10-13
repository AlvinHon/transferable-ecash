use ark_ec::pairing::Pairing;
use ark_std::rand::RngCore;
use ark_std::UniformRand;

use super::ciphertext::Ciphertext;

/// The encryption key for the encryption scheme E - ElGamal encryption.
///
/// It is a wrapper around the `bls_elgamal::EncryptKey` struct. Additionally, it
/// implements `adapt_proof` to adapt an equality proof of knowledge.
pub struct EncryptKey<E: Pairing> {
    pub(crate) inner: (
        // D1
        bls_elgamal::EncryptKey<E::G1>,
        // D2
        bls_elgamal::EncryptKey<E::G1>,
    ),
}

impl<E: Pairing> EncryptKey<E> {
    /// Encrypt a message of (m1, m2).
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
    pub fn encrypt<R: RngCore>(
        &self,
        rng: &mut R,
        m1: E::G1Affine,
        m2: E::G1Affine,
    ) -> Ciphertext<E> {
        let v = E::ScalarField::rand(rng);
        self.encrypt_with(m1, m2, v)
    }

    /// Encrypt a message of (m1, m2) with random `v`.
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
    /// let c = ek.encrypt_with(m1, m2, v);
    /// let (m1_, m2_) = dk.decrypt(&c);
    /// assert_eq!((m1, m2), (m1_, m2_));
    /// ```
    pub fn encrypt_with(
        &self,
        m1: E::G1Affine,
        m2: E::G1Affine,
        v: E::ScalarField,
    ) -> Ciphertext<E> {
        let c1 = self.inner.0.encrypt(m1, v);
        let c2 = self.inner.1.encrypt(m2, v);
        (c1, c2).into()
    }

    /// Randomize a ciphertext.
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
    /// let c = ek.rerandomize(rng, &c);
    /// let (m1_, m2_) = dk.decrypt(&c);
    /// assert_eq!((m1, m2), (m1_, m2_));
    /// ```
    pub fn rerandomize<R: RngCore>(&self, rng: &mut R, c: &Ciphertext<E>) -> Ciphertext<E> {
        let v = E::ScalarField::rand(rng);
        let (c1, c2) = c.into();
        let c1 = self.inner.0.rerandomize(c1, v);
        let c2 = self.inner.1.rerandomize(c2, v);
        (c1, c2).into()
    }

    /// Verify a ciphertext. It is equivalent to encrypting the message and comparing the ciphertexts.
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
    /// let c = ek.encrypt_with(m1, m2, v);
    /// assert!(ek.verify(m1, m2, &c, v));
    /// ```
    pub fn verify(
        &self,
        m1: E::G1Affine,
        m2: E::G1Affine,
        c: &Ciphertext<E>,
        v: E::ScalarField,
    ) -> bool {
        self.encrypt_with(m1, m2, v) == *c
    }

    pub fn adapt_proof(&self) {
        todo!()
    }
}
