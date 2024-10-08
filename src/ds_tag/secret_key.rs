use ark_ec::pairing::Pairing;
use std::ops::Mul;

use crate::params::Params;

use super::{
    message::Message,
    serial_number::{SerialNumber, SerialNumberProof},
    tag::{Tag, TagProof},
};

#[derive(Clone, PartialEq, Eq)]
pub struct SecretKey<E: Pairing> {
    pub(crate) sk: E::ScalarField,
}

impl<E: Pairing> SecretKey<E> {
    /// The serial-number generation function, on input a secret key and
    /// a nonce, outputs a serial-number component and a message which is
    /// signed by the bank using a signature scheme.
    pub fn init_serial_number(
        &self,
        params: &Params<E>,
        n: E::ScalarField,
    ) -> (SerialNumber<E>, (Message<E>, Message<E>)) {
        // m = g1^n, n = g2^(n+sk)
        let sn_0 = params.g1.mul(n).into();
        let sn_1 = params.g2.mul(self.sk + n).into();
        // M = g1^n, N = g^n
        let message_0 = Message {
            m: params.g1.mul(n).into(),
            n: params.g.mul(n).into(),
        };
        // M' = g1^sk, N' = g^sk
        let message_1 = Message {
            m: params.g1.mul(self.sk).into(),
            n: params.g.mul(self.sk).into(),
        };
        (SerialNumber { m: sn_0, n: sn_1 }, (message_0, message_1))
    }

    /// The serial-number generation function, on input a secret key and
    /// a nonce, outputs a serial-number component and a proof of well-formedness.
    pub fn generate_serial_number(
        &self,
        params: &Params<E>,
        n: E::ScalarField,
    ) -> (SerialNumber<E>, SerialNumberProof<E>) {
        // m = g1^n, n = g2^(n+sk)
        let sn_0 = params.g1.mul(n).into();
        let sn_1 = params.g2.mul(self.sk + n).into();
        (
            SerialNumber { m: sn_0, n: sn_1 },
            SerialNumberProof {
                // sn_pf = g^n
                sn_pf: params.g.mul(n).into(),
            },
        )
    }

    /// The double-spending tag function, takes as input a secret key,
    /// a nonce and a serial number, and outputs a double-spending tag
    /// and a tag proof.
    pub fn generate_tag(
        &self,
        params: &Params<E>,
        n: E::ScalarField,
        sn: &SerialNumber<E>,
    ) -> (Tag<E>, TagProof<E>) {
        // A = m^sk + h1^n, B = n^sk + h2^n (i.e. A=g1^(n sk) + h1^n, B=g2^(n+sk)*sk + h2^n)
        let tag_0 = sn.m.mul(self.sk) + params.h1.mul(n);
        let tag_1 = sn.n.mul(self.sk) + params.h2.mul(n);
        (
            Tag {
                a: tag_0.into(),
                b: tag_1.into(),
            },
            TagProof {
                t_pf: params.g.mul(n).into(),
            },
        )
    }
}
