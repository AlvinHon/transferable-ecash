use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use std::ops::Mul;

use super::signature::Signature;

pub struct SigningKey<E: Pairing> {
    pub(crate) xy: Vec<(E::ScalarField, E::ScalarField)>,
}

impl<E: Pairing> SigningKey<E> {
    /// Signs a message using the one-time linearly homomorphic structure-preserving signature.
    ///
    /// A deterministic algorithm that takes the signing key sk and the message m, and outputs a signature.
    ///
    /// # Example
    ///
    /// ```rust
    /// use ark_std::{test_rng, UniformRand};
    /// use ark_ec::pairing::Pairing;
    /// use transferable_ecash::lhsps;
    ///
    /// type E = ark_bls12_381::Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    ///
    /// let rng = &mut test_rng();
    /// let (sk, pk) = lhsps::setup::<E, _>(rng, 5);
    /// let m: Vec<G1> = (0..5).map(|_| G1::rand(rng)).collect();
    /// let sig = sk.sign(&m).unwrap();
    /// assert!(pk.verify(&m, &sig));
    /// ```
    pub fn sign(&self, m: &[E::G1Affine]) -> Result<Signature<E>, ()> {
        if self.xy.len() != m.len() {
            return Err(());
        }
        // z = Π m^xi, r = Π m^yi
        let (z, r) = m
            .iter()
            .zip(&self.xy)
            .map(|(m, (x, y))| (m.mul(x), m.mul(y)))
            .fold((E::G1Affine::zero(), E::G1Affine::zero()), |acc, m| {
                ((acc.0 + m.0).into(), (acc.1 + m.1).into())
            });

        Ok(Signature { z, r })
    }
}
