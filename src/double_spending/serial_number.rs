use ark_ec::pairing::Pairing;

#[derive(Clone, PartialEq, Eq)]
pub struct SerialNumber<E: Pairing> {
    pub(crate) m: E::G1Affine,
    pub(crate) n: E::G1Affine,
}

pub struct SerialNumberProof<E: Pairing> {
    pub(crate) sn_pf: E::G2Affine,
}
