use ark_ec::pairing::Pairing;
use groth_sahai::prover::CProof;

pub struct Ciphertext<E: Pairing> {
    pub(crate) c: Vec<E::G1Affine>,
    pub(crate) cpf_b: CProof<E>,
    pub(crate) cpf_ps: Vec<CProof<E>>,
    pub(crate) cpf_v: CProof<E>,
    pub(crate) cpf_fgh: Vec<CProof<E>>,
    pub(crate) cpf_w: CProof<E>,
}
