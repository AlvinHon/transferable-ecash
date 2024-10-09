use ark_ec::pairing::Pairing;

use super::ciphertext::Ciphertext;

pub struct DecryptKey<E: Pairing> {
    pub alpha: Vec<E::ScalarField>,
}

impl<E: Pairing> DecryptKey<E> {
    /// Decrypts a ciphertext.
    ///
    /// A deterministic decryption algorithm which takes a ciphertext, and
    /// outputs either a plaintext or an error
    pub fn decrypt(&self, c: &Ciphertext<E>) -> Result<Vec<E::G1Affine>, ()> {
        todo!()
    }
}
