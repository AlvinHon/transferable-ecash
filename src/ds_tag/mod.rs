//! This module implements the double-spending tag scheme from Appendix B.1 of
//! `Transferable E-cash: A Cleaner Model and the First Practical Instantiation`.
//!
//! There are some differences from the original paper, where there are `**` marked
//! on the comments.

use ark_ec::pairing::Pairing;
use ark_ec::AffineRepr;
use ark_std::rand::RngCore;
use ark_std::{UniformRand, Zero};
use std::ops::{Mul, Neg};

#[derive(Clone)]
pub struct Params<E: Pairing> {
    pub(crate) g1: E::G1Affine,
    pub(crate) g: E::G2Affine,
    pub(crate) g2: E::G1Affine,
    pub(crate) h1: E::G1Affine,
    pub(crate) h2: E::G1Affine,
}

impl<E: Pairing> Params<E> {
    pub fn rand<R: RngCore>(rng: &mut R) -> Self {
        Self {
            g1: E::G1Affine::rand(rng),
            g: E::G2Affine::rand(rng),
            g2: E::G1Affine::rand(rng),
            h1: E::G1Affine::rand(rng),
            h2: E::G1Affine::rand(rng),
        }
    }
}

impl<E: Pairing> Params<E> {
    pub fn key_gen<R: RngCore>(&self, rng: &mut R) -> (SecretKey<E>, PublicKey<E>) {
        let sk = E::ScalarField::rand(rng);
        let pk = self.g.mul(sk).into();
        (SecretKey { sk }, PublicKey { pk })
    }
}

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

pub struct Message<E: Pairing> {
    pub(crate) m: E::G1Affine,
    pub(crate) n: E::G2Affine,
}

#[derive(Clone, PartialEq, Eq)]
pub struct SerialNumber<E: Pairing> {
    pub(crate) m: E::G1Affine,
    pub(crate) n: E::G1Affine,
}

pub struct SerialNumberProof<E: Pairing> {
    pub(crate) sn_pf: E::G2Affine,
}

pub struct Tag<E: Pairing> {
    pub(crate) a: E::G1Affine,
    pub(crate) b: E::G1Affine,
}

pub struct TagProof<E: Pairing> {
    pub(crate) t_pf: E::G2Affine,
}

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

#[cfg(test)]
mod tests {
    use ark_ec::pairing::Pairing;
    use ark_std::test_rng;
    use ark_std::{One, UniformRand};

    use crate::ds_tag::{detect, Params};

    use super::{PublicKey, Searcher};

    type E = ark_bls12_381::Bls12_381;
    type Fr = <E as Pairing>::ScalarField;

    struct SearcherImpl {
        pks: Vec<PublicKey<E>>,
    }

    impl Searcher<E> for SearcherImpl {
        fn search<F>(&self, f: F) -> Option<PublicKey<E>>
        where
            F: Fn(&PublicKey<E>) -> bool,
        {
            self.pks.iter().find_map(|pk| f(pk).then_some(pk.clone()))
        }
    }

    #[test]
    fn test_detect() {
        let rng = &mut test_rng();
        let params = Params::<E>::rand(rng);
        let (sk, pk) = params.key_gen(rng);
        let searcher = SearcherImpl {
            pks: vec![pk.clone()],
        };

        let n = Fr::rand(rng);

        let (sn, sn_pf) = sk.generate_serial_number(&params, n);
        assert!(pk.verify_serial_number(&params, &sn, &sn_pf));
        let (tag, tag_pf) = sk.generate_tag(&params, n, &sn);

        // double-spending
        let (sn1, sn1_pf) = sk.generate_serial_number(&params, n + Fr::one());
        assert!(pk.verify_serial_number(&params, &sn1, &sn1_pf));
        let (tag1, tag1_pf) = sk.generate_tag(&params, n + Fr::one(), &sn1);

        let (double_spender, proof) = detect(
            &searcher, &params, &sn, &sn1, &tag, &tag_pf, &tag1, &tag1_pf,
        )
        .unwrap();

        assert!(double_spender == pk);
        assert!(double_spender.verify_guilt(&params, &proof));
    }

    #[test]
    fn test_verifiability() {
        use ark_std::test_rng;
        use ark_std::{One, UniformRand};

        type E = ark_bls12_381::Bls12_381;
        type Fr = <E as Pairing>::ScalarField;

        let rng = &mut test_rng();
        let params = Params::<E>::rand(rng);
        let (sk, pk) = params.key_gen(rng);
        let (sk2, pk2) = params.key_gen(rng);
        let n = Fr::rand(rng);

        let (sn, (msg0, msg1)) = sk.init_serial_number(&params, n);
        assert!(pk.verify_first_serial_number(&params, &sn, &(msg0, msg1)));

        let (another_sn, another_pf) = sk2.generate_serial_number(&params, Fr::rand(rng));
        assert!(pk2.verify_serial_number(&params, &another_sn, &another_pf));

        // generate tag from serial number generated by another key.
        let (tagx, tagx_pf) = sk.generate_tag(&params, n, &another_sn);
        assert!(pk.verify_tag(&params, &sn, &another_sn, &tagx, &tagx_pf));

        let (sn1, sn1_pf) = sk.generate_serial_number(&params, n + Fr::one());
        assert!(pk.verify_serial_number(&params, &sn1, &sn1_pf));

        // generate tag from serial number generated by another key.
        let (tagx, tagx_pf) = sk.generate_tag(&params, n + Fr::one(), &another_sn);
        assert!(pk.verify_tag(&params, &sn1, &another_sn, &tagx, &tagx_pf));
    }

    #[test]
    fn test_2_show_extractability() {
        use ark_std::test_rng;
        use ark_std::{One, UniformRand};

        type E = ark_bls12_381::Bls12_381;
        type Fr = <E as Pairing>::ScalarField;

        struct SearcherImpl {
            pks: Vec<PublicKey<E>>,
        }

        impl Searcher<E> for SearcherImpl {
            fn search<F>(&self, f: F) -> Option<PublicKey<E>>
            where
                F: Fn(&PublicKey<E>) -> bool,
            {
                self.pks.iter().find_map(|pk| f(pk).then_some(pk.clone()))
            }
        }

        let rng = &mut test_rng();
        let params = Params::<E>::rand(rng);
        let (sk, pk) = params.key_gen(rng);
        let (sk1, pk1) = (sk.clone(), pk.clone());
        let (sk2, pk2) = (sk.clone(), pk.clone());
        let searcher = SearcherImpl {
            pks: vec![pk.clone(), pk1.clone(), pk2.clone()],
        };

        let n = Fr::rand(rng);

        // T.SVfy(pk, sn, n) = 1
        let (sn, sn_pf) = sk.generate_serial_number(&params, n);
        assert!(pk.verify_serial_number(&params, &sn, &sn_pf));

        // T.SVfy(pk1, sn1, sn-pf) = T.Sfy(pk2, sn2, sn2-pf)
        let (sn1, sn1_pf) = sk1.generate_serial_number(&params, n + Fr::one()); // different nonce hence different sn
        assert!(pk1.verify_serial_number(&params, &sn1, &sn1_pf));
        let (sn2, sn2_pf) = sk2.generate_serial_number(&params, n + Fr::one() + Fr::one()); // different nonce hence different sn
        assert!(pk2.verify_serial_number(&params, &sn2, &sn2_pf));

        // T.TVfy(pk1, sn, sn1, tag1, tag1-pf) = T.TVfy(pk2, sn, sn2, tag2, tag2-pf) = 1
        let (tag1, tag1_pf) = sk1.generate_tag(&params, n, &sn1); // <- use the same nonce but different sn
        assert!(pk1.verify_tag(&params, &sn, &sn1, &tag1, &tag1_pf));
        let (tag2, tag2_pf) = sk2.generate_tag(&params, n, &sn2); // <- use the same nonce but different sn
        assert!(pk2.verify_tag(&params, &sn, &sn2, &tag2, &tag2_pf));

        assert!(sn1 != sn2);

        let (double_spender, proof) = detect(
            &searcher, &params, &sn1, &sn2, &tag1, &tag1_pf, &tag2, &tag2_pf,
        )
        .unwrap();

        assert!(double_spender == pk);
        assert!(double_spender.verify_guilt(&params, &proof));
    }
}
