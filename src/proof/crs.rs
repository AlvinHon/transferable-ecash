//! Defines the struct `CRS` (Common Reference String).

use ark_ec::pairing::Pairing;
use ark_std::{rand::RngCore, UniformRand};
use gs_ppe::CommitmentKeys;
use ndarray::{arr2, Array2};

/// Common Reference String (CRS) for GS proof.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CRS<E: Pairing> {
    pub g1_gen: E::G1Affine,
    pub g2_gen: E::G2Affine,
    pub cks: CommitmentKeys<E>,
}

impl<E: Pairing> CRS<E> {
    pub fn rand<R: RngCore>(rng: &mut R) -> Self {
        let g1_gen = E::G1Affine::rand(rng);
        let g2_gen = E::G2Affine::rand(rng);
        let cks = CommitmentKeys::<E>::setup(rng, g1_gen, g2_gen);

        CRS {
            g1_gen,
            g2_gen,
            cks,
        }
    }

    pub(crate) fn u(&self) -> Array2<E::G1> {
        arr2(&[
            [self.cks.u.0 .0.into(), self.cks.u.1 .0.into()],
            [self.cks.u.0 .1.into(), self.cks.u.1 .1.into()],
        ])
    }

    pub(crate) fn u1(&self) -> Array2<E::G1> {
        arr2(&[[self.cks.u.0 .0.into(), self.cks.u.0 .1.into()]])
    }

    pub(crate) fn u2(&self) -> Array2<E::G1> {
        arr2(&[[self.cks.u.1 .0.into(), self.cks.u.1 .1.into()]])
    }

    pub(crate) fn v(&self) -> Array2<E::G2> {
        arr2(&[
            [self.cks.v.0 .0.into(), self.cks.v.1 .0.into()],
            [self.cks.v.0 .1.into(), self.cks.v.1 .1.into()],
        ])
    }

    pub(crate) fn v1(&self) -> Array2<E::G2> {
        arr2(&[[self.cks.v.0 .0.into(), self.cks.v.0 .1.into()]])
    }

    pub(crate) fn v2(&self) -> Array2<E::G2> {
        arr2(&[[self.cks.v.1 .0.into(), self.cks.v.1 .1.into()]])
    }
}

impl<E: Pairing> AsRef<CommitmentKeys<E>> for CRS<E> {
    fn as_ref(&self) -> &CommitmentKeys<E> {
        &self.cks
    }
}
