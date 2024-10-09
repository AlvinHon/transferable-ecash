use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::AffineRepr;
use ark_std::rand::{Rng, RngCore};
use ark_std::{One, UniformRand, Zero};
use groth_sahai::prover::{CProof, Provable};
use groth_sahai::statement::PPE;
use groth_sahai::verifier::Verifiable;
use groth_sahai::{Matrix, CRS};
use std::ops::{Mul, Neg};

use crate::lhsps;

use super::ciphertext::Ciphertext;

pub struct EncryptKey<E: Pairing> {
    pub(crate) f: E::G1Affine,
    pub(crate) g: E::G1Affine, // = crs.g1_gen

    pub(crate) h_alpha: Vec<E::G1Affine>,
    pub(crate) crs: CRS<E>,
    pub(crate) sigs: (
        lhsps::signature::Signature<E>,
        lhsps::signature::Signature<E>,
    ),
}

impl<E: Pairing> EncryptKey<E> {
    /// Encrypt a message.
    ///
    /// A randomized encryption algorithm which takes as input,
    /// - a message `m`,
    /// - some randomness `v` (same length as `m`),
    ///
    /// and outputs a ciphertext.
    pub fn encrypt<R: Rng>(
        &self,
        rng: &mut R,
        m: &[E::G1Affine],
        v: &[E::G1Affine],
    ) -> Ciphertext<E> {
        let phi = E::ScalarField::rand(rng);
        // c = [c0, c1, ..., cn+1]
        //   = [f^phi, g^phi, m1^phi + h1^phi, m2^phi + h2^phi, ..., mn^phi + hn^phi]
        let mut c = vec![self.f.mul(phi).into(), self.g.mul(phi).into()];
        c.extend(
            m.iter()
                .zip(v.iter())
                .map(|(mi, vi)| (mi.mul(phi) + vi.mul(phi)).into()),
        );

        // generate gs-proof of e(g^b, g2) + e(g, g2^-b) = 0
        // TODO: trivial because b = 1. optimization point here.
        let b = E::ScalarField::one();
        // e(A, Y) + e(X, B) = e(g, g~^-b) + e(g^b, g~) = 0
        let cpf_b = gs_proof(
            rng,
            &self.crs,
            self.crs.g1_gen,
            self.crs.g2_gen.mul(b.neg()).into(),
            self.g.mul(b).into(),
            self.crs.g2_gen,
        );
        // generate gs-proof of e(ps_i, g2) + e(ci, g2^-b) = 0
        let cpf_ps: Vec<CProof<E>> = c
            .iter()
            .map(|c_i| {
                // ps_i = c_i^b
                let ps_i = c_i.mul(b).into();
                // e(A, Y) + e(X, B) = e(c_i, g~^-b) + e(g^b, g~) = 0
                gs_proof::<E, _>(
                    rng,
                    &self.crs,
                    *c_i,
                    self.crs.g2_gen.mul(b.neg()).into(),
                    ps_i,
                    self.crs.g2_gen,
                )
            })
            .collect();
        // v = [c_0^b, c_1^b, g^(1-b), c_1^(1-b), ..., c_n+1^(1-b)]
        let mut v = vec![c[0].mul(b).into(), c[1].mul(b).into()];
        v.push(self.g.mul(b.neg()).into());
        v.extend(
            c.iter()
                .skip(1)
                .map(|c_i| c_i.mul(b.neg()).into())
                .collect::<Vec<_>>(),
        );
        // generate lhsps signature on v
        // TODO

        // generate proof of validity of lhsps signature on v
        // TODO cps_v

        // generate proof of (f^b, g^b, h_1^b, ..., h_n^b)
        // TODO cps_fgh

        // w = (f^b, g^b, 1, h_1^(1-b), ..., h_n^(1-b))
        // TODO

        // generate lhsps signature on w
        // TODO

        // generate proof of validity of lhsps signature on w
        // TODO cps_w

        // Output ciphertext c = (ci for i in 1..n, cps_b, cps_ps, cps_v, cps_fgh, cps_w)
        todo!()
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
    pub fn adpt_proof(&self) {
        todo!()
    }
}

/// Create GS proof for pairing product equation: e(A, Y) + e(X, B) = 0.
/// This function is used by encryption function in EncryptKey.
pub(crate) fn gs_proof<E: Pairing, R: RngCore>(
    rng: &mut R,
    crs: &CRS<E>,
    a: E::G1Affine,
    y: E::G2Affine,
    x: E::G1Affine,
    b: E::G2Affine,
) -> CProof<E> {
    // Apply:
    //      Π e(A_i, Y_i) + Π e(X_i, B_i) + ΠΠ e(X_i, Y_j)^gamma_ij = t
    // We have:
    //  n = 1, m = 1,
    //  A = [a], B = [0, b], X = [0, x], Y = [y],
    //  gamma = 0
    let xvars = vec![E::G1Affine::zero(), x];
    let yvars = vec![y];
    let a_consts = vec![a];
    let b_consts = vec![E::G2Affine::zero(), b];
    let gamma: Matrix<E::ScalarField> =
        vec![vec![E::ScalarField::zero()], vec![E::ScalarField::zero()]];
    let target: PairingOutput<E> = PairingOutput::<E>::zero();
    let equ: PPE<E> = PPE::<E> {
        a_consts,
        b_consts,
        gamma,
        target,
    };
    let proof: CProof<E> = equ.commit_and_prove(&xvars, &yvars, &crs, rng);
    assert!(equ.verify(&proof, &crs));
    proof
}
