use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_std::rand::Rng;
use ark_std::{One, UniformRand};
use groth_sahai::prover::CProof;
use groth_sahai::CRS;
use std::ops::{Mul, Neg};

use crate::lhsps;
use crate::proof::create_proof_ayxb;

use super::ciphertext::Ciphertext;

pub struct EncryptKey<E: Pairing> {
    pub(crate) f: E::G1Affine,
    pub(crate) g: E::G1Affine,

    pub(crate) h: Vec<E::G1Affine>,
    pub(crate) crs: CRS<E>,

    pub(crate) lhsps_sig_v1: lhsps::signature::Signature<E>,
    pub(crate) lhsps_sig_v2: lhsps::signature::Signature<E>,
    pub(crate) lhsps_vk: lhsps::verifying_key::VerifyKey<E>,
}

impl<E: Pairing> EncryptKey<E> {
    /// Encrypt a message.
    ///
    /// A randomized encryption algorithm which takes as input,
    /// - a message `m`,
    /// - some randomness `rng`,
    ///
    /// and outputs a ciphertext.
    pub fn encrypt<R: Rng>(&self, rng: &mut R, m: &[E::G1Affine]) -> Ciphertext<E> {
        let phi = E::ScalarField::rand(rng);
        // c = [c0, c1, ..., cn+1]
        //   = [f^phi, g^phi, m1^phi + h1^phi, m2^phi + h2^phi, ..., mn^phi + hn^phi]
        let mut c = vec![self.f.mul(phi).into(), self.g.mul(phi).into()];
        c.extend(
            m.iter()
                .zip(self.h.iter())
                .map(|(mi, hi)| (*mi + hi.mul(phi)).into()),
        );

        // generate gs-proof of e(g^b, g2) + e(g, g2^-b) = 0
        // TODO: trivial because b = 1. optimization point here.
        let b = E::ScalarField::one();
        // e(A, Y) + e(X, B) = e(g, g~^-b) + e(g^b, g~) = 0
        let cpf_b = create_proof_ayxb(
            rng,
            &self.crs,
            self.g,
            self.crs.g2_gen.mul(b.neg()).into(),
            self.g.mul(b).into(),
            self.crs.g2_gen,
        )
        .unwrap();
        // generate gs-proof of e(ps_i, g2) + e(ci, g2^-b) = 0
        let cpf_ps: Vec<CProof<E>> = c
            .iter()
            .skip(1)
            .map(|c_i| {
                // ps_i = c_i^b
                let ps_i = c_i.mul(b).into();
                // e(A, Y) + e(X, B) = e(c_i, g~^-b) + e(g^b, g~) = 0
                create_proof_ayxb::<E, _>(
                    rng,
                    &self.crs,
                    *c_i,
                    self.crs.g2_gen.mul(b.neg()).into(),
                    ps_i,
                    self.crs.g2_gen,
                )
                .unwrap()
            })
            .collect();
        // v = [c_0^b, c_1^b, g^(1-b), c_2^(1-b), ..., c_n+1^(1-b)]
        //   = [c_0, c_1, 1, 1, ..., 1]
        let mut v = vec![
            c[0].mul(b).into(),
            c[1].mul(b).into(),
            self.g.mul(E::ScalarField::one() + b.neg()).into(),
        ];
        v.extend(
            c.iter()
                .skip(3)
                .map(|c_i| c_i.mul(E::ScalarField::one() + b.neg()).into())
                .collect::<Vec<_>>(),
        );
        // generate lhsps signature on v = v1^phi + v2 ^ 0 = v1^phi, hence only lhsps_sig_v1 is needed
        let sig_with_w = vec![(phi, self.lhsps_sig_v1)];
        let sigv = self.lhsps_vk.sign_derive(&sig_with_w).unwrap();

        // generate proof of validity of lhsps signature on v
        let cpf_v = self
            .lhsps_vk
            .generate_proof(rng, &self.crs, &v, &sigv)
            .unwrap();

        // generate proof of (f^b, g^b, h_1^b, ..., h_n^b)
        let mut fgh = vec![self.f.mul(b).into(), self.g.mul(b).into()];
        fgh.extend(
            self.h
                .iter()
                .map(|h_i| h_i.mul(b).into())
                .collect::<Vec<_>>(),
        );
        // generate gs-proof of:
        // 1. e(f^b, g2) + e(f, g2^-b) = 0
        // 2. e(g^b, g2) + e(g, g2^-b) = 0 (i.e. cpf_b)
        // 3. e(h_i^b, g2) + e(h_i, g2^-b) = 0
        let cpf_fgh: Vec<CProof<E>> = fgh
            .iter()
            .map(|fgh_i| {
                // e(A, Y) + e(X, B) = 0
                create_proof_ayxb::<E, _>(
                    rng,
                    &self.crs,
                    *fgh_i,
                    self.crs.g2_gen.mul(b.neg()).into(),
                    fgh_i.mul(b).into(),
                    self.crs.g2_gen,
                )
                .unwrap()
            })
            .collect();

        // w = (f^b, g^b, 1, h_1^(1-b), ..., h_n^(1-b))
        //   = (f, g, 1, 1, ..., 1)
        let mut w = vec![
            self.f.mul(b).into(),
            self.g.mul(b).into(),
            E::G1Affine::zero(),
        ];
        w.extend(
            self.h
                .iter()
                .map(|h_i| h_i.mul(E::ScalarField::one() + b.neg()).into())
                .collect::<Vec<_>>(),
        );

        // generate lhsps signature on w = v1^b + v2^0 = v1^b, hence only lhsps_sig_v1 is needed
        let sig_with_w = vec![(b, self.lhsps_sig_v1)];
        let sigw = self.lhsps_vk.sign_derive(&sig_with_w).unwrap();

        // generate proof of validity of lhsps signature on w
        let cpf_w = self
            .lhsps_vk
            .generate_proof(rng, &self.crs, &w, &sigw)
            .unwrap();

        // Output ciphertext c = (ci for i in 1..n, cpf_b, cpf_ps, cpf_v, cpf_fgh, cpf_w)
        Ciphertext {
            c,
            cpf_b,
            cpf_ps,
            cpf_v,
            cpf_fgh,
            cpf_w,
        }
    }

    /// Re-randomize a ciphertext.
    ///
    /// A randomized algorithm which mutates the input ciphertext `c` with some randomness `v`.
    pub fn rerandomize(&self, c: &mut Ciphertext<E>, v: &[E::G1Affine]) {
        todo!()
    }

    /// Verify a ciphertext.
    ///
    /// A deterministic algorithm which takes as input,
    /// - a message `m`,
    /// - some randomness `v` (same length as `m`),
    /// - a ciphertext `c`m
    ///
    /// and outputs a bit.
    pub fn verify(&self, m: &[E::G1Affine], v: &[E::G1Affine], c: &Ciphertext<E>) -> bool {
        todo!()
    }

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
