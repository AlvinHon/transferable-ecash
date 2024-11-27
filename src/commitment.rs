//! Defines the functions of the commit-and-prove systems (as notated as `C`) in section 4.2 of the paper.

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_std::{rand::Rng, One, Zero};
use gs_ppe::{Com, CommitmentKeys, Equation, Matrix, Proof, ProofSystem, Variable};
use std::ops::Mul;

use crate::double_spending::{self as T, serial_number::SerialNumberProof};

pub(crate) struct C<E: Pairing> {
    cks: CommitmentKeys<E>,
}

impl<E: Pairing> C<E> {
    /// `C.Setup(Gr)`
    pub(crate) fn setup<R: ark_std::rand::RngCore>(rng: &mut R) -> Self {
        Self {
            cks: CommitmentKeys::rand(rng),
        }
    }

    // ... Auxiliary functions ...

    /// `C.Prv_{sn,init}`
    pub(crate) fn prove_init_serial_number<R: Rng>(
        &self,
        rng: &mut R,
        ds_params: &T::DSParams<E>,
        pk: &T::PublicKey<E>,
        sn: &T::SerialNumber<E>,
        msgs: &(T::Message<E>, T::Message<E>),
    ) -> SnInitProof<E> {
        // Constructs the GS proof of the following equations (which are used in the `verify_first_serial_number` function):
        // e(N, g) + e(g2^-1, M1) + e(g2^-1, pk) == 1
        // e(M1, g) == e(g1, M1)
        // e(M2, g) == e(g1, M2)

        let T::DSParams { g, g1, g2, .. } = ds_params;
        let pk = pk.pk;

        // 1.
        let var_n = Variable::new(rng, sn.n);
        let var_m1 = Variable::new(rng, msgs.0.n);
        let var_pk = Variable::new(rng, pk);
        let neg_g2 = g2.mul(-E::ScalarField::one()).into();
        let ProofSystem {
            c: c1,
            d: d1,
            proof: proof1,
            ..
        } = gs_ppe::setup(
            rng,
            &self.cks,
            &[(neg_g2, var_m1), (neg_g2, var_pk)],
            &[(var_n, *g)],
            &Matrix::new(&[[E::ScalarField::zero(), E::ScalarField::zero()]]), // dim = 1x2
        );

        // 2.
        let var_m1_y = Variable::new(rng, msgs.0.n);
        let var_m1_x = Variable::new(rng, msgs.0.m);
        let neg_g = g.mul(-E::ScalarField::one()).into();
        let ProofSystem {
            c: c2,
            d: d2,
            proof: proof2,
            ..
        } = gs_ppe::setup(
            rng,
            &self.cks,
            &[(*g1, var_m1_y)],
            &[(var_m1_x, neg_g)],
            &Matrix::new(&[[E::ScalarField::zero()]]), // dim = 1x1
        );

        // 3.
        let var_m2_y = Variable::new(rng, msgs.1.n);
        let var_m2_x = Variable::new(rng, msgs.1.m);
        let neg_g = g.mul(-E::ScalarField::one()).into();
        let ProofSystem {
            c: c3,
            d: d3,
            proof: proof3,
            ..
        } = gs_ppe::setup(
            rng,
            &self.cks,
            &[(*g1, var_m2_y)],
            &[(var_m2_x, neg_g)],
            &Matrix::new(&[[E::ScalarField::zero()]]), // dim = 1x1
        );

        SnInitProof {
            c1,
            d1,
            proof1,
            c2,
            d2,
            proof2,
            c3,
            d3,
            proof3,
        }
    }

    /// `C.Verify_{sn,init}`
    pub(crate) fn verify_init_serial_number(
        &self,
        ds_params: &T::DSParams<E>,
        proof: &SnInitProof<E>,
    ) -> bool {
        let T::DSParams { g, g1, g2, .. } = ds_params;

        // 1.
        let neg_g2 = g2.mul(-E::ScalarField::one()).into();
        let eq1 = Equation::<E>::new(
            vec![neg_g2, neg_g2],
            vec![*g],
            Matrix::new(&[[E::ScalarField::zero(), E::ScalarField::zero()]]),
            PairingOutput::zero(),
        );

        if !eq1.verify(&self.cks, &proof.c1, &proof.d1, &proof.proof1) {
            return false;
        }

        // 2.
        let neg_g = g.mul(-E::ScalarField::one()).into();
        let eq2 = Equation::<E>::new(
            vec![*g1],
            vec![neg_g],
            Matrix::new(&[[E::ScalarField::zero()]]),
            PairingOutput::zero(),
        );

        if !eq2.verify(&self.cks, &proof.c2, &proof.d2, &proof.proof2) {
            return false;
        }

        // 3.
        let eq3 = Equation::<E>::new(
            vec![*g1],
            vec![neg_g],
            Matrix::new(&[[E::ScalarField::zero()]]),
            PairingOutput::zero(),
        );

        eq3.verify(&self.cks, &proof.c3, &proof.d3, &proof.proof3)
    }

    /// `C.Prv_sn`
    pub(crate) fn prove_serial_number<R: Rng>(
        &self,
        rng: &mut R,
        ds_params: &T::DSParams<E>,
        pk: &T::PublicKey<E>,
        sn: &T::SerialNumber<E>,
        sn_pf: &SerialNumberProof<E>,
    ) -> SnProof<E> {
        // Constructs the GS proof of the following equations (which are used in the `verify_first_serial_number` function):
        // e(N, g) + e(g2^-1, sn-pf) + e(g2^-1, pk) == 1
        // e(M, g) + e(g1^-1, sn-pf) == 1

        let T::DSParams { g, g1, g2, .. } = ds_params;
        let pk = pk.pk;

        // 1.
        let var_n = Variable::new(rng, sn.n);
        let var_sn_pf = Variable::new(rng, sn_pf.sn_pf);
        let var_pk = Variable::new(rng, pk);
        let neg_g2 = g2.mul(-E::ScalarField::one()).into();
        let ProofSystem {
            c: c1,
            d: d1,
            proof: proof1,
            ..
        } = gs_ppe::setup(
            rng,
            &self.cks,
            &[(neg_g2, var_sn_pf), (neg_g2, var_pk)],
            &[(var_n, *g)],
            &Matrix::new(&[[E::ScalarField::zero(), E::ScalarField::zero()]]), // dim = 1x2
        );

        // 2.
        let var_m = Variable::new(rng, sn.m);
        let var_sn_pf = Variable::new(rng, sn_pf.sn_pf);
        let neg_g1 = g1.mul(-E::ScalarField::one()).into();
        let ProofSystem {
            c: c2,
            d: d2,
            proof: proof2,
            ..
        } = gs_ppe::setup(
            rng,
            &self.cks,
            &[(neg_g1, var_sn_pf)],
            &[(var_m, *g)],
            &Matrix::new(&[[E::ScalarField::zero()]]), // dim = 1x1
        );

        SnProof {
            c1,
            d1,
            proof1,
            c2,
            d2,
            proof2,
        }
    }

    /// `C.Verify_sn`
    pub(crate) fn verify_serial_number(
        &self,
        ds_params: &T::DSParams<E>,
        proof: &SnProof<E>,
    ) -> bool {
        let T::DSParams { g, g1, g2, .. } = ds_params;

        // 1.
        let neg_g2 = g2.mul(-E::ScalarField::one()).into();
        let eq1 = Equation::<E>::new(
            vec![neg_g2, neg_g2],
            vec![*g],
            Matrix::new(&[[E::ScalarField::zero(), E::ScalarField::zero()]]),
            PairingOutput::zero(),
        );

        if !eq1.verify(&self.cks, &proof.c1, &proof.d1, &proof.proof1) {
            return false;
        }

        // 2.
        let neg_g1 = g1.mul(-E::ScalarField::one()).into();
        let eq2 = Equation::<E>::new(
            vec![neg_g1],
            vec![*g],
            Matrix::new(&[[E::ScalarField::zero()]]),
            PairingOutput::zero(),
        );

        eq2.verify(&self.cks, &proof.c2, &proof.d2, &proof.proof2)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct SnInitProof<E: Pairing> {
    // e(N, g) + e(g2^-1, M1) + e(g2^-1, pk) == 1
    pub(crate) c1: Vec<Com<<E as Pairing>::G1>>,
    pub(crate) d1: Vec<Com<<E as Pairing>::G2>>,
    pub(crate) proof1: Proof<E>,

    // e(M1, g) == e(g1, M1)
    pub(crate) c2: Vec<Com<<E as Pairing>::G1>>,
    pub(crate) d2: Vec<Com<<E as Pairing>::G2>>,
    pub(crate) proof2: Proof<E>,

    // e(M2, g) == e(g1, M2)
    pub(crate) c3: Vec<Com<<E as Pairing>::G1>>,
    pub(crate) d3: Vec<Com<<E as Pairing>::G2>>,
    pub(crate) proof3: Proof<E>,
}

#[derive(Clone, Debug)]
pub(crate) struct SnProof<E: Pairing> {
    // e(N, g) + e(g2^-1, sn-pf) + e(g2^-1, pk) == 1
    pub(crate) c1: Vec<Com<<E as Pairing>::G1>>,
    pub(crate) d1: Vec<Com<<E as Pairing>::G2>>,
    pub(crate) proof1: Proof<E>,

    // e(M, g) + e(g1^-1, sn-pf) == 1
    pub(crate) c2: Vec<Com<<E as Pairing>::G1>>,
    pub(crate) d2: Vec<Com<<E as Pairing>::G2>>,
    pub(crate) proof2: Proof<E>,
}

#[cfg(test)]
mod tests {
    use ark_ec::pairing::Pairing;
    use ark_std::test_rng;
    use ark_std::UniformRand;

    use super::*;
    use crate::double_spending::{key_gen, DSParams};

    type E = ark_bls12_381::Bls12_381;
    type Fr = <E as Pairing>::ScalarField;

    #[test]
    fn test_c_prv_sn_init() {
        let rng = &mut test_rng();

        let ds_params = DSParams::<E>::rand(rng);
        let (sk, pk) = key_gen(rng, &ds_params);

        let n = Fr::rand(rng);

        let (sn, msgs) = sk.init_serial_number(&ds_params, n);
        assert!(pk.verify_first_serial_number(&ds_params, &sn, &msgs));

        let c = C::<E>::setup(rng);

        let proof = c.prove_init_serial_number(rng, &ds_params, &pk, &sn, &msgs);
        assert!(c.verify_init_serial_number(&ds_params, &proof));
    }

    #[test]
    fn test_c_prv_sn() {
        let rng = &mut test_rng();

        let ds_params = DSParams::<E>::rand(rng);
        let (sk, pk) = key_gen(rng, &ds_params);

        let n = Fr::rand(rng);

        let (sn, sn_pf) = sk.generate_serial_number(&ds_params, n);
        assert!(pk.verify_serial_number(&ds_params, &sn, &sn_pf));

        let c = C::<E>::setup(rng);

        let proof = c.prove_serial_number(rng, &ds_params, &pk, &sn, &sn_pf);
        assert!(c.verify_serial_number(&ds_params, &proof));
    }
}
