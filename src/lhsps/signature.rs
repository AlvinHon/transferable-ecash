use ark_ec::pairing::Pairing;

pub struct Signature<E: Pairing> {
    pub(crate) z: E::G1Affine,
    pub(crate) r: E::G1Affine,
}
