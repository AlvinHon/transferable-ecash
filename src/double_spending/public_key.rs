use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_std::Zero;
use std::ops::Neg;

use crate::params::Params;

use super::{
    detect::DetectionProof,
    message::Message,
    serial_number::{SerialNumber, SerialNumberProof},
    tag::{Tag, TagProof},
};

#[derive(Clone, PartialEq, Eq)]
pub struct PublicKey<E: Pairing> {
    pub(crate) pk: E::G2Affine,
}

impl<E: Pairing> PublicKey<E> {
    pub fn from(pk: E::G2Affine) -> Self {
        Self { pk }
    }

    /// The incrimination-proof verification function.
    pub fn verify_guilt(&self, params: &Params<E>, proof: &DetectionProof<E>) -> bool {
        // ** this is different from the version that I read. I guess there is a missing part on the paper. **
        // ** Original: e(ax, g) == e(mx, pk) **
        // e(ax, g) == e(mx, pk) + e(hx, tx)
        let lhs = E::pairing(proof.ax, params.g);
        let rhs_0 = E::pairing(proof.mx, self.pk);
        lhs == rhs_0 + E::pairing(params.h1, proof.tx)
            || lhs == rhs_0 + E::pairing(params.h2, proof.tx)
    }

    /// On input a public key, a serial number and a message, checks their consistency.
    pub fn verify_first_serial_number(
        &self,
        params: &Params<E>,
        sn: &SerialNumber<E>,
        msgs: &(Message<E>, Message<E>),
    ) -> bool {
        // e(M, g) + e(g1^-1, M1) == 1
        let eq = E::pairing(sn.m, params.g) + E::pairing(params.g1.into_group().neg(), msgs.0.n);
        if !eq.is_zero() {
            return false;
        }
        // ** this is different from the version that I read. I guess there is typo on the paper. **
        // ** Original: e(M, g) + e(g2^-1, M2) + e(g2^-1, pk) == 1 **
        // e(N, g) + e(g2^-1, M1) + e(g2^-1, pk) == 1
        let eq = E::pairing(sn.n, params.g)
            + E::pairing(params.g2.into_group().neg(), msgs.0.n)
            + E::pairing(params.g2.into_group().neg(), self.pk);
        if !eq.is_zero() {
            return false;
        }
        // M2 == pk
        if msgs.1.n != self.pk {
            return false;
        }
        // e(M1, g) == e(g1, M1)
        let eq = E::pairing(msgs.0.m, params.g) - E::pairing(params.g1, msgs.0.n);
        if !eq.is_zero() {
            return false;
        }
        // e(M2, g) == e(g1, M2)
        let eq = E::pairing(msgs.1.m, params.g) - E::pairing(params.g1, msgs.1.n);
        if !eq.is_zero() {
            return false;
        }

        true
    }

    /// On input a public key and a serial number, checks their consistency.
    pub fn verify_serial_number(
        &self,
        params: &Params<E>,
        sn: &SerialNumber<E>,
        sn_pf: &SerialNumberProof<E>,
    ) -> bool {
        // e(M, g) + e(g1^-1, sn-pf) == 1
        let eq = E::pairing(sn.m, params.g) + E::pairing(params.g1.into_group().neg(), sn_pf.sn_pf);
        if !eq.is_zero() {
            return false;
        }
        // e(N, g) + e(g2^-1, sn-pf) + e(g2^-1, pk) == 1
        let eq = E::pairing(sn.n, params.g)
            + E::pairing(params.g2.into_group().neg(), sn_pf.sn_pf)
            + E::pairing(params.g2.into_group().neg(), self.pk);
        if !eq.is_zero() {
            return false;
        }

        true
    }

    /// On input a public key, two serial numbers, a double-spending tag, and a proof,
    /// checks consistency of the tag w.r.t the key and the serial numbers.
    pub fn verify_tag(
        &self,
        params: &Params<E>,
        sn: &SerialNumber<E>,
        sn_d: &SerialNumber<E>,
        tag: &Tag<E>,
        tag_pf: &TagProof<E>,
    ) -> bool {
        // e(M, g) + e(g1^-1, tag-pf) == 1
        let eq = E::pairing(sn.m, params.g) + E::pairing(params.g1.into_group().neg(), tag_pf.t_pf);
        if !eq.is_zero() {
            return false;
        }
        // e(A, g^-1) + e(M_d, pk) + e(h1, tag-pf) == 1
        let eq = E::pairing(tag.a, params.g.into_group().neg())
            + E::pairing(sn_d.m, self.pk)
            + E::pairing(params.h1, tag_pf.t_pf);
        if !eq.is_zero() {
            return false;
        }
        // e(B, g^-1) + e(N_d, pk) + e(h2, tag-pf) == 1
        let eq = E::pairing(tag.b, params.g.into_group().neg())
            + E::pairing(sn_d.n, self.pk)
            + E::pairing(params.h2, tag_pf.t_pf);
        if !eq.is_zero() {
            return false;
        }

        true
    }
}
