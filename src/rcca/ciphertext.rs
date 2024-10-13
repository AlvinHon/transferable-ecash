use ark_ec::{pairing::Pairing, AffineRepr};

use groth_sahai::prover::CProof;

use crate::proof::check_proof_ayxb;

use super::encrypt_key::EncryptKey;

pub struct Ciphertext<E: Pairing> {
    pub(crate) c: Vec<E::G1Affine>,
    pub(crate) cpf_b: CProof<E>,
    pub(crate) cpf_ps: Vec<CProof<E>>,
    pub(crate) cpf_v: CProof<E>,
    pub(crate) cpf_fgh: Vec<CProof<E>>,
    pub(crate) cpf_w: CProof<E>,
}

impl<E: Pairing> Ciphertext<E> {
    /// Check all proofs of the ciphertext.
    pub fn check_proofs(&self, enc_key: &EncryptKey<E>) -> Result<(), ()> {
        if self.cpf_ps.len() != self.c.len() - 1 {
            return Err(());
        }

        // check all proofs
        let crs = &enc_key.crs;
        // cfp_b is proof of e(A, Y) + e(X, B) = e(g, g~^-b) + e(g^b, g~) = 0
        if !check_proof_ayxb(crs, &self.cpf_b, enc_key.g, crs.g2_gen) {
            return Err(());
        }
        // cfp_ps is proof of e(A, Y) + e(X, B) = e(c_i, g~^-b) + e(g^b, g~) = 0
        if self
            .c
            .iter()
            .skip(1)
            .zip(self.cpf_ps.iter())
            .any(|(ci, cpf)| !check_proof_ayxb(crs, cpf, *ci, crs.g2_gen))
        {
            return Err(());
        }
        // cpf_v is proof for message v = [c_0, c_1, 1, ..., 1]
        let mut v = vec![self.c[0], self.c[1]];
        v.extend(vec![E::G1Affine::zero(); self.c.len() - 2]);
        if !enc_key.lhsps_vk.check_proof(crs, &self.cpf_v, &v) {
            return Err(());
        }
        // cpf_fgh is proof for message fgh = (f, g, h_1, ..., h_n)
        let mut fgh = vec![enc_key.f, enc_key.g];
        fgh.extend(enc_key.h.iter());
        if fgh
            .iter()
            .zip(self.cpf_fgh.iter())
            .any(|(fgh_i, cpf)| !check_proof_ayxb(crs, cpf, *fgh_i, crs.g2_gen))
        {
            return Err(());
        }
        // cpf_w is proof for message w = [f, g, 1, 1, ..., 1]
        let mut w = vec![enc_key.f, enc_key.g];
        w.extend(vec![E::G1Affine::zero(); self.c.len() - 2]);
        if !enc_key.lhsps_vk.check_proof(crs, &self.cpf_w, &w) {
            return Err(());
        }

        Ok(())
    }
}
