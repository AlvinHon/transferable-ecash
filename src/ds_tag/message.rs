use ark_ec::pairing::Pairing;

pub struct Message<E: Pairing> {
    pub(crate) m: E::G1Affine,
    pub(crate) n: E::G2Affine,
}
