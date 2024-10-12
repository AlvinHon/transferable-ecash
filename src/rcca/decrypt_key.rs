use ark_ec::{pairing::Pairing, AffineRepr};
use std::ops::{Mul, Neg};

use crate::proof::check_proof_ayxb;

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
        if c.cpf_ps.len() != c.c.len() - 1 {
            return Err(());
        }

        // check all proofs
        let crs = &self.enc_key.crs;
        // cfp_b is proof of e(A, Y) + e(X, B) = e(g, g~^-b) + e(g^b, g~) = 0
        if !check_proof_ayxb(crs, &c.cpf_b, self.enc_key.g, crs.g2_gen) {
            return Err(());
        }
        // cfp_ps is proof of e(A, Y) + e(X, B) = e(c_i, g~^-b) + e(g^b, g~) = 0
        if c.c
            .iter()
            .skip(1)
            .zip(c.cpf_ps.iter())
            .any(|(ci, cpf)| !check_proof_ayxb(crs, cpf, *ci, crs.g2_gen))
        {
            return Err(());
        }
        // cpf_v is proof for message v = [c_0, c_1, 1, ..., 1]
        let mut v = vec![c.c[0], c.c[1]];
        v.extend(vec![E::G1Affine::zero(); c.c.len() - 2]);
        if !self.enc_key.lhsps_vk.check_proof(crs, &c.cpf_v, &v) {
            return Err(());
        }
        // cpf_fgh is proof for message fgh = (f, g, h_1, ..., h_n)
        let mut fgh = vec![self.enc_key.f, self.enc_key.g];
        fgh.extend(self.enc_key.h.iter());
        if fgh
            .iter()
            .zip(c.cpf_fgh.iter())
            .any(|(fgh_i, cpf)| !check_proof_ayxb(crs, cpf, *fgh_i, crs.g2_gen))
        {
            return Err(());
        }
        // cpf_w is proof for message w = [f, g, 1, 1, ..., 1]
        let mut w = vec![self.enc_key.f, self.enc_key.g];
        w.extend(vec![E::G1Affine::zero(); c.c.len() - 2]);
        if !self.enc_key.lhsps_vk.check_proof(crs, &c.cpf_w, &w) {
            return Err(());
        }

        // compute M_i = c_i+1 / c_1^alpha_i
        let mut m = Vec::new();
        for i in 0..c.c.len() - 2 {
            m.push((c.c[i + 2] + c.c[1].mul(self.alpha[i]).neg()).into());
        }

        Ok(m)
    }
}
