use ark_ec::pairing::Pairing;

#[derive(Clone, PartialEq, Eq)]
pub struct Ciphertext<E: Pairing> {
    pub c0: E::G1Affine,
    pub c1: E::G1Affine,
    pub c2: E::G1Affine,
}

impl<'a, E: Pairing>
    From<(
        bls_elgamal::Ciphertext<E::G1>,
        bls_elgamal::Ciphertext<E::G1>,
    )> for Ciphertext<E>
{
    fn from(
        c: (
            bls_elgamal::Ciphertext<E::G1>,
            bls_elgamal::Ciphertext<E::G1>,
        ),
    ) -> Ciphertext<E> {
        Ciphertext {
            c0: c.0 .0,
            c1: c.0 .1,
            c2: c.1 .1,
        }
    }
}

impl<'a, E: Pairing> From<&'a Ciphertext<E>>
    for (
        bls_elgamal::Ciphertext<E::G1>,
        bls_elgamal::Ciphertext<E::G1>,
    )
{
    fn from(
        c: &'a Ciphertext<E>,
    ) -> (
        bls_elgamal::Ciphertext<E::G1>,
        bls_elgamal::Ciphertext<E::G1>,
    ) {
        (
            bls_elgamal::Ciphertext(c.c0, c.c1),
            bls_elgamal::Ciphertext(c.c0, c.c2),
        )
    }
}
