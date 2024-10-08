use ark_ec::pairing::Pairing;

pub struct Tag<E: Pairing> {
    pub(crate) a: E::G1Affine,
    pub(crate) b: E::G1Affine,
}

pub struct TagProof<E: Pairing> {
    pub(crate) t_pf: E::G2Affine,
}
