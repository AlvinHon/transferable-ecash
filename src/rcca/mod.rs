//! This module implements the Replayable-CCA encryption scheme from Appendix B.2 of
//! of `Transferable E-cash: A Cleaner Model and the First Practical Instantiation`.

use std::vec;

use ark_ec::{pairing::Pairing, AffineRepr};
use ark_std::{ops::Mul, rand::RngCore, UniformRand};
use decrypt_key::DecryptKey;
use encrypt_key::EncryptKey;
use groth_sahai::{AbstractCrs, CRS};

use crate::lhsps;

pub mod ciphertext;
pub mod decrypt_key;
pub mod encrypt_key;

pub fn key_gen<E: Pairing, R: RngCore>(rng: &mut R, n: usize) -> (DecryptKey<E>, EncryptKey<E>) {
    let crs = CRS::<E>::generate_crs(rng);
    let crs_cloned = // TODO derive Clone for CRS
        CRS::<E> {
            u: crs.u.clone(),
            v: crs.v.clone(),
            g1_gen: crs.g1_gen,
            g2_gen: crs.g2_gen,
            gt_gen: crs.gt_gen,
        };

    let f = E::G1Affine::rand(rng);
    let g = E::G1Affine::rand(rng);

    // alpha = [alpha1, alpha2, ..., alphan]
    let alpha = (0..n)
        .map(|_| E::ScalarField::rand(rng))
        .collect::<Vec<_>>();
    // hi = g^alpha_i for i in 1..n
    let h = alpha
        .iter()
        .map(|alpha_i| g.mul(alpha_i).into())
        .collect::<Vec<_>>();
    // **
    // v1 = [f,g,1,1,...,1]
    let mut v1 = vec![f, g];
    v1.extend(vec![E::G1Affine::zero(); n + 1]);
    // v2 = [1,1,1,h1,h2,...,hn]
    let mut v2 = vec![E::G1Affine::zero(); 3];
    v2.extend_from_slice(&h);

    // Notice that the LHSPS signing key tk will never be published by the key
    // generation algorithm, it will only be used in the security proofs.
    let (tk, lhsps_vk) = lhsps::setup::<E, _>(rng, n + 3);
    let lhsps_sig_v1 = tk.sign(&v1).unwrap();
    let lhsps_sig_v2 = tk.sign(&v2).unwrap();

    (
        DecryptKey {
            enc_key: EncryptKey {
                f,
                g,
                h: h.clone(),
                crs: crs_cloned,
                lhsps_sig_v1,
                lhsps_sig_v2,
                lhsps_vk: lhsps_vk.clone(),
            },

            alpha,
        },
        EncryptKey {
            f,
            g,
            h,
            crs,
            lhsps_sig_v1,
            lhsps_sig_v2,
            lhsps_vk,
        },
    )
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_std::UniformRand;

    use crate::rcca::key_gen;

    type E = Bls12_381;
    type G1 = <E as Pairing>::G1Affine;

    #[test]
    fn debug_test() {
        let rng = &mut ark_std::test_rng();
        let (sk, pk) = key_gen::<E, _>(rng, 5);
        let m = (0..5).map(|_| G1::rand(rng)).collect::<Vec<_>>();
        let c = pk.encrypt(rng, &m);
        let m_d = sk.decrypt(&c).unwrap();
        assert_eq!(m, m_d);
    }
}
