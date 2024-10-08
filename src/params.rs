use ark_ec::pairing::Pairing;
use ark_std::rand::RngCore;
use ark_std::UniformRand;

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
