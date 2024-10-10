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
    /// use ark_std::{test_rng, UniformRand};
    /// use ark_ec::pairing::Pairing;
    /// use transferable_ecash::lhsps;
    /// use std::ops::Mul;
    ///
    /// type E = ark_bls12_381::Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    /// type Fr = <E as Pairing>::ScalarField;
    ///
    /// let rng = &mut test_rng();
    /// let (sk, pk) = lhsps::setup::<E, _>(rng, 5);
    /// let mut m: Vec<G1> = (0..5).map(|_| G1::rand(rng)).collect();
    /// let sig = sk.sign(&m).unwrap();
    /// let w = Fr::rand(rng);
    /// let sig_with_w = vec![(w, sig)];
    /// let sig_d = pk.sign_derive(&sig_with_w).unwrap();
    /// let m_d = m.iter().map(|mi| mi.mul(w).into()).collect::<Vec<_>>();
    /// assert!(pk.verify(&m_d, &sig_d));
    /// ```
    pub fn sign_derive(
        &self,
        sig_with_w: &[(E::ScalarField, Signature<E>)],
    ) -> Result<Signature<E>, ()> {
        if sig_with_w.is_empty() {
            return Err(());
        }
        // cannot derive signature on message with Mi' = Mi^0 = 1. (verify must fail)
        if sig_with_w.iter().any(|(w, _)| w.is_zero()) {
            return Err(());
        }

        // z = Π z^w, r = Π r^w
        let (z, r) = sig_with_w
            .iter()
            .map(|(w, sigi)| (sigi.z.mul(*w), sigi.r.mul(*w)))
            .fold((E::G1Affine::zero(), E::G1Affine::zero()), |acc, p| {
                ((acc.0 + p.0).into(), (acc.1 + p.1).into())
            });

        Ok(Signature { z, r })
    }

    /// Verifies a signature using the one-time linearly homomorphic structure-preserving signature.
    ///
    /// A deterministic algorithm taking the verification key pk, the message vector m and a signature.
    /// It outputs true if the signature is valid, false otherwise.
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
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use groth_sahai::{AbstractCrs, CRS};
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
    ///
    /// let crs = CRS::<E>::generate_crs(rng);
    /// pk.generate_proof(rng, &crs, &m, &sig)
    ///     .expect("proof should be valid");
    /// ```
    pub fn generate_proof<R: RngCore>(
        &self,
        rng: &mut R,
        crs: &CRS<E>,
        m: &[E::G1Affine],
        sig: &Signature<E>,
    ) -> Result<CProof<E>, ()> {
        let target = m
            .iter()
            .zip(&self.pk)
            .map(|(m, pk)| E::pairing(*m, *pk))
            .fold(PairingOutput::zero(), |acc, m| acc + m);
        gs_proof_xbxb_t(rng, &crs, sig.z, self.gz, sig.r, self.gr, target)
    }
}

#[cfg(test)]
mod tests {
    use crate::lhsps::setup;
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_std::UniformRand;
    use groth_sahai::{AbstractCrs, CRS};
    use std::ops::Mul;

    type E = Bls12_381;
    type G1 = <E as Pairing>::G1Affine;
    type Fr = <E as Pairing>::ScalarField;

    #[test]
    fn test_derived_signature() {
        let rng = &mut ark_std::test_rng();
        let crs = CRS::<E>::generate_crs(rng);

        let (sk, pk) = setup::<E, _>(rng, 5);
        let m = (0..5).map(|_| G1::rand(rng)).collect::<Vec<_>>();
        let sig = sk.sign(&m).unwrap();
        pk.generate_proof(rng, &crs, &m, &sig)
            .expect("proof on original message should be valid");

        // derive signature from original signature
        let w = Fr::rand(rng);
        let sig_d = pk.sign_derive(&[(w, sig)]).unwrap();
        let m_d = m.iter().map(|mi| mi.mul(w).into()).collect::<Vec<_>>();
        assert!(pk.verify(&m_d, &sig_d));

        // proof generated from derived signature should be valid on derived message m_d
        pk.generate_proof(rng, &crs, &m_d, &sig_d)
            .expect("proof should be valid");
    }

    #[test]
    fn test_two_derived_signatures() {
        let rng = &mut ark_std::test_rng();
        let crs = CRS::<E>::generate_crs(rng);

        let (sk, pk) = setup::<E, _>(rng, 5);
        let m1 = (0..5).map(|_| G1::rand(rng)).collect::<Vec<_>>();
        let sig1 = sk.sign(&m1).unwrap();
        pk.generate_proof(rng, &crs, &m1, &sig1)
            .expect("proof on original message 1 should be valid");

        let m2 = (0..5).map(|_| G1::rand(rng)).collect::<Vec<_>>();
        let sig2 = sk.sign(&m2).unwrap();
        pk.generate_proof(rng, &crs, &m2, &sig2)
            .expect("proof on original message 2 should be valid");

        // derive signature from original signature
        let w1 = Fr::rand(rng);
        let w2 = Fr::rand(rng);
        let sig_d = pk.sign_derive(&[(w1, sig1), (w2, sig2)]).unwrap();
        let m1_d = m1.iter().map(|mi| mi.mul(w1).into()).collect::<Vec<G1>>();
        let m2_d = m2.iter().map(|mi| mi.mul(w2).into()).collect::<Vec<G1>>();
        let m1_m2_d = m1_d
            .iter()
            .zip(&m2_d)
            .map(|(m1i, m2i)| (*m1i + *m2i).into())
            .collect::<Vec<_>>();
        assert!(pk.verify(&m1_m2_d, &sig_d));

        // proof generated from derived signature should be valid on derived message m1_d + m2_d
        pk.generate_proof(rng, &crs, &m1_m2_d, &sig_d)
            .expect("proof should be valid");
    }
}
