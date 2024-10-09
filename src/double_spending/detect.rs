use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use std::ops::Neg;

use super::{
    public_key::PublicKey,
    serial_number::SerialNumber,
    tag::{Tag, TagProof},
};
use crate::params::Params;

pub fn detect<E: Pairing, S: Searcher<E>>(
    searcher: &S,
    params: &Params<E>,
    sn0: &SerialNumber<E>,
    sn1: &SerialNumber<E>,
    tag0: &Tag<E>,
    tag0_pf: &TagProof<E>,
    tag1: &Tag<E>,
    tag1_pf: &TagProof<E>,
) -> Option<(PublicKey<E>, DetectionProof<E>)> {
    let a = (tag0.a + tag1.a.into_group().neg()).into();
    let b = (tag0.b + tag1.b.into_group().neg()).into();
    let m = (sn0.m + sn1.m.into_group().neg()).into();
    let n = (sn0.n + sn1.n.into_group().neg()).into();
    let tx = (tag0_pf.t_pf + tag1_pf.t_pf.into_group().neg()).into();

    let (ax, mx, hx) = if a == E::G1Affine::zero() {
        (b, n, params.h2)
    } else {
        (a, m, params.h1)
    };

    searcher
        .search(|pk| {
            // ** this is different from the version that I read. I guess there is a missing part on the paper. **
            // ** Original: e(A, g) == e(M, pk) **
            // e(A, g) == e(M, pk) + e(hx, t)
            E::pairing(ax, params.g) == E::pairing(mx, pk.pk) + E::pairing(hx, tx)
        })
        .map(|pk| (pk, DetectionProof { ax, mx, tx }))
}

pub trait Searcher<E: Pairing> {
    fn search<F>(&self, f: F) -> Option<PublicKey<E>>
    where
        F: Fn(&PublicKey<E>) -> bool;
}

pub struct DetectionProof<E: Pairing> {
    pub(crate) ax: E::G1Affine,
    pub(crate) mx: E::G1Affine,
    // ** this is different from the version that I read. I guess there is a missing part on the paper. **
    // ** Originally, this variable is not included. **
    pub(crate) tx: E::G2Affine,
}
