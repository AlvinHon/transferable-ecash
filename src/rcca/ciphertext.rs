use ark_ec::pairing::Pairing;
use groth_sahai::prover::CProof;

pub struct Ciphertext<E: Pairing> {
    pub(crate) c: Vec<E::G1Affine>,
    pub(crate) cps_b: CProof<E>,
    pub(crate) cps_ps: Vec<CProof<E>>,
    pub(crate) cps_v: CProof<E>,
    pub(crate) cps_fgh: CProof<E>,
    pub(crate) cps_w: CProof<E>,
}
