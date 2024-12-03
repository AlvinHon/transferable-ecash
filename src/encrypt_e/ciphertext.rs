use ark_ec::pairing::Pairing;

#[derive(Clone, PartialEq, Eq)]
pub struct Ciphertext<E: Pairing> {
    pub c0: E::G1Affine,
    pub c1: E::G1Affine,
    pub c2: E::G1Affine,
}
