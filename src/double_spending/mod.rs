//! This module implements the double-spending tag scheme from Appendix B.1 of
//! `Transferable E-cash: A Cleaner Model and the First Practical Instantiation`.
//!
//! There are some differences from the original paper, where there are `**` marked
//! on the comments.

pub mod detect;
pub mod message;
pub use message::Message;
pub mod params;
pub use params::DSParams;
pub mod public_key;
pub use public_key::PublicKey;
pub mod secret_key;
pub use secret_key::SecretKey;
pub mod serial_number;
pub use serial_number::SerialNumber;
pub mod tag;

use ark_ec::pairing::Pairing;
use ark_std::rand::RngCore;
use ark_std::UniformRand;
use std::ops::Mul;

pub fn key_gen<E: Pairing, R: RngCore>(
    rng: &mut R,
    params: &DSParams<E>,
) -> (SecretKey<E>, PublicKey<E>) {
    let sk = E::ScalarField::rand(rng);
    let pk = params.g.mul(sk).into();
    (SecretKey { sk }, PublicKey { pk })
}

#[cfg(test)]
mod tests {
    use ark_ec::pairing::Pairing;
    use ark_std::test_rng;
    use ark_std::{One, UniformRand};

    use crate::double_spending::detect::detect;
    use crate::double_spending::{key_gen, DSParams};

    use super::detect::Searcher;
    use super::PublicKey;

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
        let params = DSParams::<E>::rand(rng);
        let (sk, pk) = key_gen(rng, &params);
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
        let params = DSParams::<E>::rand(rng);
        let (sk, pk) = key_gen(rng, &params);
        let (sk2, pk2) = key_gen(rng, &params);
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
        let rng = &mut test_rng();
        let params = DSParams::<E>::rand(rng);
        let (sk, pk) = key_gen(rng, &params);
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
