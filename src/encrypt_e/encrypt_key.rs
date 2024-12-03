use ark_ec::pairing::Pairing;
use ark_std::{rand::Rng, One};
use std::ops::Mul;

use crate::proof::{check_mse_proof_ayxb, create_mse_proof_ayxb, MSEProof, CRS};

use super::ciphertext::Ciphertext;

/// The encryption key for the encryption scheme E - ElGamal encryption.
///
/// It is a wrapper around the `bls_elgamal::EncryptKey` struct. Additionally, it
/// implements `adapt_proof` to adapt an equality proof of knowledge.
pub struct EncryptKey<E: Pairing> {
    /// The group generator.
    pub(crate) g: E::G1Affine,
    pub(crate) y1: E::G1Affine, // g^dk1
    pub(crate) y2: E::G1Affine, // g^dk2
}

impl<E: Pairing> EncryptKey<E> {
    /// Encrypt a message of (m1, m2) with randomness `v`.
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
    ///
    /// let c = ek.encrypt(m1, m2, v);
    /// let (m1_, m2_) = dk.decrypt(&c);
    /// assert_eq!((m1, m2), (m1_, m2_));
    /// ```
    pub fn encrypt(&self, m1: E::G1Affine, m2: E::G1Affine, v: E::ScalarField) -> Ciphertext<E> {
        let c0 = self.g.mul(v).into();
        let c1 = (self.y1.mul(v) + m1).into();
        let c2 = (self.y2.mul(v) + m2).into();
        Ciphertext { c0, c1, c2 }
    }

    /// Create a proof of knowledge of the messages `m1` and `m2`, and the randomness `v` used in the encryption.
    ///
    /// # Example
    ///
    /// ```rust
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use transferable_ecash::{encrypt_e, proof};
    ///
    /// type E = ark_bls12_381::Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    /// type Fr = <E as Pairing>::ScalarField;
    ///
    /// let rng = &mut test_rng();
    /// let crs = proof::CRS::<E>::rand(rng);
    /// let (dk, ek) = encrypt_e::key_gen::<E, _>(rng);
    /// let (m1, m2) = (G1::rand(rng), G1::rand(rng));
    /// let v = Fr::rand(rng);
    ///
    /// let c = ek.encrypt(m1, m2, v);
    /// let proofs = ek.prove(rng, &crs, m1, m2, v);
    /// assert!(ek.verify(&crs, &c, &proofs));
    /// ```
    pub fn prove<R: Rng>(
        &self,
        rng: &mut R,
        crs: &CRS<E>,
        m1: E::G1Affine,
        m2: E::G1Affine,
        v: E::ScalarField,
    ) -> (MSEProof<E>, MSEProof<E>, MSEProof<E>) {
        // prove knowledge of v s.t. c0 = g^v
        let proof1 = create_mse_proof_ayxb(rng, crs, vec![self.g.into()], vec![v], vec![], vec![]);

        // prove knowledge of v and m1 s.t. c1 = m1 + y1^v
        let proof2 = create_mse_proof_ayxb(
            rng,
            crs,
            vec![self.y1.into()],
            vec![v],
            vec![m1.into()],
            vec![E::ScalarField::one()],
        );

        // prove knowledge of v and m2 s.t. c2 = m2 + y2^v
        let proof3 = create_mse_proof_ayxb(
            rng,
            crs,
            vec![self.y2.into()],
            vec![v],
            vec![m2.into()],
            vec![E::ScalarField::one()],
        );

        (proof1, proof2, proof3)
    }

    /// Randomize a ciphertext with randomness `v`.
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
    /// let v_prime = Fr::rand(rng);
    ///
    /// let c = ek.encrypt(m1, m2, v);
    /// let c = ek.rerandomize(&c, v_prime);
    /// let (m1_, m2_) = dk.decrypt(&c);
    /// assert_eq!((m1, m2), (m1_, m2_));
    /// ```
    pub fn rerandomize(&self, c: &Ciphertext<E>, v: E::ScalarField) -> Ciphertext<E> {
        let c0 = (c.c0 + self.g.mul(v)).into();
        let c1 = (c.c1 + self.y1.mul(v)).into();
        let c2 = (c.c2 + self.y2.mul(v)).into();
        Ciphertext { c0, c1, c2 }
    }

    /// Verify proof that proves the knowledge of messages `m1` and `m2`, and the randomness `v`
    /// used in the encryption associated with the ciphertext.
    ///
    /// # Example
    ///
    /// ```rust
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use transferable_ecash::{encrypt_e, proof};
    ///
    /// type E = ark_bls12_381::Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    /// type Fr = <E as Pairing>::ScalarField;
    ///
    /// let rng = &mut test_rng();
    /// let crs = proof::CRS::<E>::rand(rng);
    /// let (dk, ek) = encrypt_e::key_gen::<E, _>(rng);
    /// let (m1, m2) = (G1::rand(rng), G1::rand(rng));
    /// let v = Fr::rand(rng);
    ///
    /// let c = ek.encrypt(m1, m2, v);
    /// let proofs = ek.prove(rng, &crs, m1, m2, v);
    /// assert!(ek.verify(&crs, &c, &proofs));
    /// ```
    pub fn verify(
        &self,
        crs: &CRS<E>,
        c: &Ciphertext<E>,
        proofs: &(MSEProof<E>, MSEProof<E>, MSEProof<E>),
    ) -> bool {
        check_mse_proof_ayxb(crs, vec![self.g.into()], vec![], c.c0.into(), &proofs.0)
            && check_mse_proof_ayxb(
                crs,
                vec![self.y1.into()],
                vec![E::ScalarField::one()],
                c.c1.into(),
                &proofs.1,
            )
            && check_mse_proof_ayxb(
                crs,
                vec![self.y2.into()],
                vec![E::ScalarField::one()],
                c.c2.into(),
                &proofs.2,
            )
    }

    /// Adapt a proof to a rerandomization (with randomness `v` used in ciphertext re-randomization)
    /// and outputs equality proofs.
    ///
    /// # Example
    ///
    /// ```rust
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::{test_rng, UniformRand};
    /// use transferable_ecash::{encrypt_e, proof};
    ///
    /// type E = ark_bls12_381::Bls12_381;
    /// type G1 = <E as Pairing>::G1Affine;
    /// type Fr = <E as Pairing>::ScalarField;
    ///
    /// let rng = &mut test_rng();
    /// let crs = proof::CRS::<E>::rand(rng);
    /// let (dk, ek) = encrypt_e::key_gen::<E, _>(rng);
    /// let (m1, m2) = (G1::rand(rng), G1::rand(rng));
    /// let v = Fr::rand(rng);
    ///
    /// let c = ek.encrypt(m1, m2, v);
    /// let proofs = ek.prove(rng, &crs, m1, m2, v);
    ///
    /// let v_prime = Fr::rand(rng);
    /// let c_prime = ek.rerandomize(&c, v_prime);
    /// let proofs_prime = ek.adapt_proof(rng, &crs, proofs, v_prime);
    /// assert!(ek.verify(&crs, &c_prime, &proofs_prime));
    /// ```
    pub fn adapt_proof<R: Rng>(
        &self,
        rng: &mut R,
        crs: &CRS<E>,
        proofs: (MSEProof<E>, MSEProof<E>, MSEProof<E>), // included the commitments
        v: E::ScalarField,
    ) -> (MSEProof<E>, MSEProof<E>, MSEProof<E>) {
        let (proof1, proof2, proof3) = proofs;
        let proof1 = proof1.adapt_proof(crs, vec![v], vec![]).randomize(
            rng,
            crs,
            vec![self.g.into()],
            vec![],
        );
        let proof2 = proof2.adapt_proof(crs, vec![v], vec![]).randomize(
            rng,
            crs,
            vec![self.y1.into()],
            vec![E::ScalarField::one()],
        );
        let proof3 = proof3.adapt_proof(crs, vec![v], vec![]).randomize(
            rng,
            crs,
            vec![self.y2.into()],
            vec![E::ScalarField::one()],
        );

        (proof1, proof2, proof3)
    }
}
