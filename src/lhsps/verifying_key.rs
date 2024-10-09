use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::AffineRepr;
use ark_std::rand::RngCore;
use ark_std::Zero;
use groth_sahai::prover::CProof;
use groth_sahai::CRS;
use std::ops::Mul;

use crate::proof::gs_proof_xbxb_t;

use super::signature::Signature;

pub struct VerifyKey<E: Pairing> {
    pub(crate) gz: E::G2Affine,
    pub(crate) gr: E::G2Affine,

    pub(crate) pk: Vec<E::G2Affine>,
}

impl<E: Pairing> VerifyKey<E> {
    /// Derives a new signature and message from the given signature and message.
    ///
    /// A deterministic algorithm taking the verification key pk and signature pairs (wi, sigi) for i in 1..l.
    /// It outputs a signature and a new message.
    ///
    /// # Example
    ///
    /// ```rust
    /// use ark_bls12_381::Bls12_381;
    /// use ark_std::test_rng;
    /// use ark_std::UniformRand;
    /// use ark_std::rand::RngCore;
    /// use ark_ec::pairing::Pairing;
    /// use transferable_ecash::lhsps;
    ///
    /// type E = Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    /// type Fr = <E as Pairing>::ScalarField;
    ///
    /// let rng = &mut test_rng();
    /// let (sk, pk) = lhsps::setup::<E, _>(rng, 5);
    /// let mut m: Vec<G1> = (0..5).map(|_| G1::rand(rng)).collect();
    /// let sig = sk.sign(&m).unwrap();
    /// let sig_with_w = vec![(Fr::rand(rng), sig)];
    /// let (m_d, sig_d) = pk.sign_derive(vec![m].as_ref(), &sig_with_w).unwrap();
    /// assert!(pk.verify(&m_d, &sig_d));
    /// ```
    pub fn sign_derive(
        &self,
        m: &[Vec<E::G1Affine>],
        sig_with_w: &[(E::ScalarField, Signature<E>)],
    ) -> Result<(Vec<E::G1Affine>, Signature<E>), ()> {
        if sig_with_w.len() != m.len() || m.iter().any(|mi| mi.len() != self.pk.len()) {
            return Err(());
        }

        let l = sig_with_w.len();
        let n = self.pk.len();

        // m = (Πm1^w1 for _ in 1..l, Πm2^w2 for _ in 1..l, ..., Πmn^wl for _ in 1..m,)
        let mut new_m = Vec::new();
        for j in 0..n {
            let mut m_sum = E::G1Affine::zero();
            for i in 0..l {
                let w = sig_with_w[i].0;
                m_sum = (m_sum + m[i][j].mul(w)).into();
            }
            new_m.push(m_sum);
        }

        // z = Π z^w, r = Π r^w
        let (z, r) = sig_with_w
            .iter()
            .map(|(w, sigi)| (sigi.z.mul(*w), sigi.r.mul(*w)))
            .fold((E::G1Affine::zero(), E::G1Affine::zero()), |acc, p| {
                ((acc.0 + p.0).into(), (acc.1 + p.1).into())
            });

        Ok((new_m, Signature { z, r }))
    }

    /// Verifies a signature using the one-time linearly homomorphic structure-preserving signature.
    ///
    /// A deterministic algorithm taking the verification key pk, the message vector m and a signature.
    /// It outputs true if the signature is valid, false otherwise.
    ///
    /// # Example
    ///
    /// ```rust
    /// use ark_bls12_381::Bls12_381;
    /// use ark_std::test_rng;
    /// use ark_std::UniformRand;
    /// use ark_std::rand::RngCore;
    /// use ark_ec::pairing::Pairing;
    /// use transferable_ecash::lhsps;
    ///
    /// type E = Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    ///
    /// let rng = &mut test_rng();
    /// let (sk, pk) = lhsps::setup::<E, _>(rng, 5);
    /// let m: Vec<G1> = (0..5).map(|_| G1::rand(rng)).collect();
    /// let sig = sk.sign(&m).unwrap();
    /// assert!(pk.verify(&m, &sig));
    /// ```
    pub fn verify(&self, m: &[E::G1Affine], sig: &Signature<E>) -> bool {
        if m.len() != self.pk.len() {
            return false;
        }

        // check (M1, ...,Mn) != (1, ... ,1)
        if m.iter().all(|mi| mi.is_zero()) {
            return false;
        }

        // e(z, gz)e(r, gr) == Π e(m, pk)
        let lhs = E::pairing(sig.r, self.gr) + E::pairing(sig.z, self.gz);
        let rhs = m
            .iter()
            .zip(&self.pk)
            .map(|(m, pk)| E::pairing(*m, *pk))
            .fold(PairingOutput::zero(), |acc, m| acc + m);
        lhs == rhs
    }

    /// Create GS proof from a signature that satisfies pairing product equation: e(z, gz)e(r, gr) == Π e(m, pk).
    ///
    /// # Example
    ///
    /// ```rust
    /// use ark_bls12_381::Bls12_381;
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::rand::RngCore;
    /// use ark_std::test_rng;
    /// use ark_std::UniformRand;
    /// use transferable_ecash::lhsps;
    /// use groth_sahai::{AbstractCrs, CRS};
    ///
    /// type E = Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    ///
    /// let rng = &mut test_rng();
    /// let (sk, pk) = lhsps::setup::<E, _>(rng, 5);
    /// let m: Vec<G1> = (0..5).map(|_| G1::rand(rng)).collect();
    /// let sig = sk.sign(&m).unwrap();
    /// assert!(pk.verify(&m, &sig));
    ///
    /// let crs = CRS::<E>::generate_crs(rng);
    /// let pf = pk.generate_proof(rng, &crs, &m, &sig);
    /// assert!(!pf.equ_proofs.is_empty());
    /// ```
    pub fn generate_proof<R: RngCore>(
        &self,
        rng: &mut R,
        crs: &CRS<E>,
        m: &[E::G1Affine],
        sig: &Signature<E>,
    ) -> CProof<E> {
        let target = m
            .iter()
            .zip(&self.pk)
            .map(|(m, pk)| E::pairing(*m, *pk))
            .fold(PairingOutput::zero(), |acc, m| acc + m);
        gs_proof_xbxb_t(rng, &crs, sig.z, self.gz, sig.r, self.gr, target)
    }
}
