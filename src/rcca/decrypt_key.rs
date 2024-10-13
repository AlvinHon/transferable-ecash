use ark_ec::pairing::Pairing;
use std::ops::{Mul, Neg};

use super::{ciphertext::Ciphertext, encrypt_key::EncryptKey};

pub struct DecryptKey<E: Pairing> {
    // used for proof verification in decryption.
    pub(crate) enc_key: EncryptKey<E>,

    // secret key
    pub(crate) alpha: Vec<E::ScalarField>,
}

impl<E: Pairing> DecryptKey<E> {
    /// Decrypts a ciphertext.
    ///
    /// A deterministic decryption algorithm which takes a ciphertext, and
    /// outputs either a plaintext or an error
    pub fn decrypt(&self, c: &Ciphertext<E>) -> Result<Vec<E::G1Affine>, ()> {
        // check all proofs
        c.check_proofs(&self.enc_key)?;

        // compute M_i = c_i+1 / c_1^alpha_i
        let mut m = Vec::new();
        for i in 0..c.c.len() - 2 {
            m.push((c.c[i + 2] + c.c[1].mul(self.alpha[i]).neg()).into());
        }

        Ok(m)
    }
}
