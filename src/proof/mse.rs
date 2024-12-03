//! This module provides functions related to GS proof for multi-scalar multiplication equation.

use ark_ec::{
    pairing::{Pairing, PairingOutput},
    CurveGroup,
};
use ark_std::{rand::Rng, One, UniformRand, Zero};
use ndarray::{arr2, Array2, Axis};
use std::ops::Mul;

use super::CRS;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MSEProof<E: Pairing> {
    c: Array2<E::G1>,
    d: Array2<E::G2>,
    pi: Array2<E::G2>,
    theta: Array2<E::G1>,
}

impl<E: Pairing> MSEProof<E> {
    pub(crate) fn randomize<R: Rng>(
        &self,
        rng: &mut R,
        crs: &CRS<E>,
        a: Vec<E::G1>,
        b: Vec<E::ScalarField>,
    ) -> Self {
        let m = self.c.dim().0;
        let n = self.d.dim().0;

        let r = Array2::from_shape_fn((m, 2), |_| E::ScalarField::rand(rng));
        let s = Array2::from_shape_fn((n, 1), |_| E::ScalarField::rand(rng));

        let new_c = randomize_com_x(&crs, &r, &self.c);
        let new_d = randomize_com_y(&crs, &s, &self.d);

        let t = Array2::from_shape_fn((1, 2), |_| E::ScalarField::rand(rng));

        let a = Array2::from_shape_vec((n, 1), a).unwrap();
        let b = Array2::from_shape_vec((m, 1), b).unwrap();

        let (pi, theta) = randomize_proof(crs, &r, &s, &t, &a, &b, &self.pi, &self.theta);

        Self {
            c: new_c,
            d: new_d,
            pi,
            theta,
        }
    }

    /// Adapt the proof to new changes in the variables. It enables re-randomization of the message being proved.
    pub(crate) fn adapt_proof(
        &self,
        crs: &CRS<E>,
        y_delta: Vec<E::ScalarField>,
        x_delta: Vec<E::G1>,
    ) -> Self {
        let mut c = self.c.clone();
        let mut d = self.d.clone();

        if !y_delta.is_empty() {
            assert!(y_delta.len() == d.dim().0);
            let y_delta = Array2::from_shape_vec((y_delta.len(), 1), y_delta).unwrap();
            d += &lz2(crs, &y_delta);
        }

        if !x_delta.is_empty() {
            assert!(x_delta.len() == c.dim().0);
            let x_delta = Array2::from_shape_vec((x_delta.len(), 1), x_delta).unwrap();
            c += &l1(&x_delta);
        }

        MSEProof {
            c,
            d,
            pi: self.pi.clone(),
            theta: self.theta.clone(),
        }
    }
}

/// Create a GS proof for the multi-scalar multiplication equation yA + bX = c on G1.
pub(crate) fn create_mse_proof_ayxb<E: Pairing, R: Rng>(
    rng: &mut R,
    crs: &CRS<E>,
    a: Vec<E::G1>,
    y: Vec<E::ScalarField>,
    x: Vec<E::G1>,
    b: Vec<E::ScalarField>,
) -> MSEProof<E> {
    let m = x.len();
    let n = y.len();
    assert!(m == b.len());
    assert!(n == a.len());

    let a = Array2::from_shape_vec((n, 1), a).unwrap();
    let y = Array2::from_shape_vec((n, 1), y).unwrap();
    let x = Array2::from_shape_vec((m, 1), x).unwrap();
    let b = Array2::from_shape_vec((m, 1), b).unwrap();

    let r = Array2::from_shape_fn((m, 2), |_| E::ScalarField::rand(rng));
    let s = Array2::from_shape_fn((n, 1), |_| E::ScalarField::rand(rng));
    let t = Array2::from_shape_fn((1, 2), |_| E::ScalarField::rand(rng));

    let c = commit_x(crs, &r, &x);
    let d = commit_y(crs, &s, &y);

    let (pi, theta) = proof(crs, &r, &s, &t, &a, &b);

    MSEProof { c, d, pi, theta }
}

/// Verify a GS proof for the multi-scalar multiplication equation yA + bX = c on G1.
pub(crate) fn check_mse_proof_ayxb<E: Pairing>(
    crs: &CRS<E>,
    a: Vec<E::G1>,
    b: Vec<E::ScalarField>,
    target: E::G1,
    proof: &MSEProof<E>,
) -> bool {
    let m = b.len();
    let n = a.len();

    let MSEProof { c, d, pi, theta } = proof;
    let a = Array2::from_shape_vec((n, 1), a).unwrap();
    let b = Array2::from_shape_vec((m, 1), b).unwrap();

    // l(a) d + c l(b) = l_t(target) + u pi + F(theta, v1)
    let lhs = matmul::<E>(&l1(&a).reversed_axes(), d)
        + matmul::<E>(&c.clone().reversed_axes(), &lz2(crs, &b));
    let rhs =
        l_t(crs, target) + matmul::<E>(&crs.u().reversed_axes(), pi) + f::<E>(theta, &crs.v1());

    lhs == rhs
}

/// Commit variable (group element) x in GS Proof.
fn commit_x<E: Pairing>(
    crs: &CRS<E>,
    r: &Array2<E::ScalarField>,
    x: &Array2<E::G1>,
) -> Array2<E::G1> {
    // c = l(x) + Ru
    l1(x) + scalar_matmul_g1::<E>(r, &crs.u())
}

fn randomize_com_x<E: Pairing>(
    crs: &CRS<E>,
    r: &Array2<E::ScalarField>, // m x 2
    c: &Array2<E::G1>,
) -> Array2<E::G1> {
    // c' = c + Ru
    c + scalar_matmul_g1::<E>(&r, &crs.u())
}

/// Commit variable (scalar) y in GS Proof.
fn commit_y<E: Pairing>(
    crs: &CRS<E>,
    s: &Array2<E::ScalarField>,
    y: &Array2<E::ScalarField>,
) -> Array2<E::G2> {
    // d = l(y) + s v1
    lz2(crs, y) + scalar_matmul_g2::<E>(s, &crs.v1())
}

fn randomize_com_y<E: Pairing>(
    crs: &CRS<E>,
    s: &Array2<E::ScalarField>, // n x 1
    d: &Array2<E::G2>,
) -> Array2<E::G2> {
    // d' = d + s v1
    d + scalar_matmul_g2::<E>(&s, &crs.v1())
}

/// Create a GS proof. Returns pi and theta.
fn proof<E: Pairing>(
    crs: &CRS<E>,
    r: &Array2<E::ScalarField>, // m x 2
    s: &Array2<E::ScalarField>, // n x 1
    t: &Array2<E::ScalarField>, // 1 x 2
    a: &Array2<E::G1>,          // n x 1
    b: &Array2<E::ScalarField>, // m x 1
) -> (Array2<E::G2>, Array2<E::G1>) {
    // phi = R^T l(b) - T^T v1
    let phi = scalar_matmul_g2::<E>(&r.clone().reversed_axes(), &lz2(crs, b))
        - scalar_matmul_g2::<E>(&t.clone().reversed_axes(), &crs.v1());

    // theta = s^T l(a) + T u
    let theta = scalar_matmul_g1::<E>(&s.clone().reversed_axes(), &l1(a))
        + scalar_matmul_g1::<E>(t, &crs.u());

    (phi, theta)
}

fn randomize_proof<E: Pairing>(
    crs: &CRS<E>,
    r: &Array2<E::ScalarField>, // m x 2
    s: &Array2<E::ScalarField>, // n x 1
    t: &Array2<E::ScalarField>, // 1 x 2
    a: &Array2<E::G1>,          // n x 1
    b: &Array2<E::ScalarField>, // m x 1
    pi: &Array2<E::G2>,
    theta: &Array2<E::G1>,
) -> (Array2<E::G2>, Array2<E::G1>) {
    // phi' = phi + R^T l(b) - T^T v1
    let phi = pi.clone() + scalar_matmul_g2::<E>(&r.clone().reversed_axes(), &lz2(crs, b))
        - scalar_matmul_g2::<E>(&t.clone().reversed_axes(), &crs.v1());

    // theta' = theta + s^T l(a) + T u
    let theta = theta.clone()
        + scalar_matmul_g1::<E>(&s.clone().reversed_axes(), &l1(a))
        + scalar_matmul_g1::<E>(t, &crs.u());

    (phi, theta)
}

/// Mapping function of Group elements, l(a) = [0 | a], dim = (m, 2)
fn l1<G: CurveGroup>(a: &Array2<G>) -> Array2<G> {
    let a = a.clone();
    let m = a.dim().0;
    let mut zeros = Array2::from_elem((m, 1), G::zero());
    zeros.append(Axis(1), a.view()).unwrap(); // dim = (m, 2)
    zeros
}

/// Mapping function of Field elements, // l(z) = z u where u = u2 + (0, p). dim = (m, 2)
fn lz1<E: Pairing>(crs: &CRS<E>, z: &Array2<E::ScalarField>) -> Array2<E::G1> {
    let m = z.dim().0;

    // u = u2 + (0, p) = [u[0, 1] | u[1, 1] + p]
    let mut u = Array2::from_elem((m, 1), (crs.cks.u.1 .0).into());
    let ps = Array2::from_elem((m, 1), crs.cks.u.1 .1 + crs.g1_gen);
    u.append(Axis(1), ps.view()).unwrap(); // dim = (m, 2)

    // l(z) = z u
    u * z
}

/// Mapping function of Field elements, // l(z) = z v where v = v2 + (0, p). dim = (m, 2)
fn lz2<E: Pairing>(crs: &CRS<E>, z: &Array2<E::ScalarField>) -> Array2<E::G2> {
    let m = z.dim().0;

    // v = v2 + (0, p) = [v[0, 1] | v[1, 1] + p]
    let mut v = Array2::from_elem((m, 1), (crs.cks.v.1 .0).into());
    let ps = Array2::from_elem((m, 1), crs.cks.v.1 .1 + crs.g2_gen);
    v.append(Axis(1), ps.view()).unwrap(); // dim = (m, 2)

    // l(z) = z u
    v * z
}

/// Mapping function of target group element G1. Returns matrix with dim = (2, 2)
fn l_t<E: Pairing>(crs: &CRS<E>, target: E::G1) -> Array2<PairingOutput<E>> {
    let x = arr2(&[[E::G1::zero(), target]]);
    let y = lz2(crs, &arr2(&[[E::ScalarField::one()]]));
    f(&x, &y)
}

/// Mapping function of paring product on group elements G1 and G2. Returns matrix with dim = (2, 2)
fn f<E: Pairing>(x: &Array2<E::G1>, y: &Array2<E::G2>) -> Array2<PairingOutput<E>> {
    // ([[x1, x2]], [[y1, y2]]) -> [[e(x1, y1), e(x1, y2)], [e(x2, y1), e(x2, y2)]]

    arr2(&[
        [
            E::pairing(x[(0, 0)], y[(0, 0)]),
            E::pairing(x[(0, 0)], y[(0, 1)]),
        ],
        [
            E::pairing(x[(0, 1)], y[(0, 0)]),
            E::pairing(x[(0, 1)], y[(0, 1)]),
        ],
    ])
}

/// Scalar matrix multiplication, left matrix is scalar and right matrix is group element (G1).
pub(crate) fn scalar_matmul_g1<E: Pairing>(
    a: &Array2<E::ScalarField>,
    b: &Array2<E::G1>,
) -> Array2<E::G1> {
    let (m, n_prime) = a.dim();
    let (m_prime, n) = b.dim();
    assert!(n_prime == m_prime);

    let mut res = Array2::from_elem((m, n), E::G1::zero());
    for i in 0..m {
        for j in 0..n {
            let mut sum = E::G1::zero();
            for k in 0..n_prime {
                sum += b[[k, j]].mul(a[[i, k]]);
            }
            res[[i, j]] = sum;
        }
    }

    res
}

/// Scalar matrix multiplication, left matrix is group element (G1) and right matrix is group element (G2).
fn matmul<E: Pairing>(a: &Array2<E::G1>, b: &Array2<E::G2>) -> Array2<PairingOutput<E>> {
    let (m, n_prime) = a.dim();
    let (m_prime, n) = b.dim();
    assert!(n_prime == m_prime);

    let mut res = Array2::from_elem((m, n), PairingOutput::zero());
    for i in 0..m {
        for j in 0..n {
            let mut sum = PairingOutput::zero();
            for k in 0..n_prime {
                sum += E::pairing(a[[i, k]], b[[k, j]]);
            }
            res[[i, j]] = sum;
        }
    }

    res
}

/// Scalar matrix multiplication, left matrix is scalar and right matrix is group element (G2).
fn scalar_matmul_g2<E: Pairing>(a: &Array2<E::ScalarField>, b: &Array2<E::G2>) -> Array2<E::G2> {
    let (m, n_prime) = a.dim();
    let (m_prime, n) = b.dim();
    assert!(n_prime == m_prime);

    let mut res = Array2::from_elem((m, n), E::G2::zero());
    for i in 0..m {
        for j in 0..n {
            let mut sum = E::G2::zero();
            for k in 0..n_prime {
                sum += b[[k, j]].mul(a[[i, k]]);
            }
            res[[i, j]] = sum;
        }
    }

    res
}

#[cfg(test)]
mod test {
    use ark_bls12_381::Bls12_381 as E;

    type G1 = <E as Pairing>::G1;
    type Fr = <E as Pairing>::ScalarField;

    use super::*;

    #[test]
    fn test_mse_proof() {
        let rng = &mut ark_std::test_rng();

        let crs = CRS::<E>::rand(rng);

        // c = m + rY
        let m = G1::rand(rng);
        let y = G1::rand(rng);
        let r = Fr::rand(rng);
        let c = m + y.mul(r);

        // GS proof for multi-scalar multiplication equation:
        // yA + bX = c
        // where y = [r], A = [Y], X = [m], b = [1]

        let proof = create_mse_proof_ayxb(rng, &crs, vec![y], vec![r], vec![m], vec![Fr::one()]);

        assert!(check_mse_proof_ayxb(
            &crs,
            vec![y],
            vec![Fr::one()],
            c,
            &proof
        ));

        // Test randomization

        let new_proof = proof.randomize(rng, &crs, vec![y], vec![Fr::one()]);
        assert!(proof.c != new_proof.c);
        assert!(proof.d != new_proof.d);
        assert!(proof.pi != new_proof.pi);
        assert!(proof.theta != new_proof.theta);

        assert!(check_mse_proof_ayxb(
            &crs,
            vec![y],
            vec![Fr::one()],
            c,
            &new_proof
        ));
    }

    #[test]
    fn test_adapt_mse_proof() {
        let rng = &mut ark_std::test_rng();

        let crs = CRS::<E>::rand(rng);

        // c = m + rY
        let m = G1::rand(rng);
        let y = G1::rand(rng);
        let r = Fr::rand(rng);
        let c = m + y.mul(r);

        // GS proof for multi-scalar multiplication equation:
        // yA + bX = c
        // where y = [r], A = [Y], X = [m], b = [1]

        let proof = create_mse_proof_ayxb(rng, &crs, vec![y], vec![r], vec![m], vec![Fr::one()]);

        assert!(check_mse_proof_ayxb(
            &crs,
            vec![y],
            vec![Fr::one()],
            c,
            &proof
        ));

        // Modify variable y
        let r_prime = Fr::rand(rng);
        let c_prime = c + y.mul(r_prime); // c' = c + r'Y = m + (r + r')Y

        let new_proof = proof.adapt_proof(&crs, vec![r_prime], vec![]);
        assert!(check_mse_proof_ayxb(
            &crs,
            vec![y],
            vec![Fr::one()],
            c_prime,
            &new_proof
        ));

        // Modify variable m
        let m_prime = G1::rand(rng);
        let c_prime = c + m_prime; // c' = m' + rY = (m + m') + rY

        let new_proof = proof.adapt_proof(&crs, vec![], vec![m_prime]);
        assert!(check_mse_proof_ayxb(
            &crs,
            vec![y],
            vec![Fr::one()],
            c_prime,
            &new_proof
        ));
    }

    #[test]
    fn test_mse_proof_ay() {
        let rng = &mut ark_std::test_rng();

        let crs = CRS::<E>::rand(rng);

        // c = rg
        let y = G1::rand(rng);
        let r = Fr::rand(rng);
        let c = y.mul(r);

        // GS proof for multi-scalar multiplication equation:
        // yA = c
        // where y = [r], A = [Y], X = [], b = []

        let proof = create_mse_proof_ayxb(rng, &crs, vec![y], vec![r], vec![], vec![]);

        assert!(check_mse_proof_ayxb(&crs, vec![y], vec![], c, &proof));

        // Test randomization

        let new_proof = proof.randomize(rng, &crs, vec![y], vec![]);
        assert!(proof.c == new_proof.c); // no commitment for x
        assert!(proof.d != new_proof.d);
        assert!(proof.pi != new_proof.pi);
        assert!(proof.theta != new_proof.theta);

        assert!(check_mse_proof_ayxb(&crs, vec![y], vec![], c, &new_proof));
    }

    #[test]
    fn test_adapt_mse_proof_ay() {
        let rng = &mut ark_std::test_rng();

        let crs = CRS::<E>::rand(rng);

        // c = rY
        let y = G1::rand(rng);
        let r = Fr::rand(rng);
        let c = y.mul(r);

        // GS proof for multi-scalar multiplication equation:
        // yA = c
        // where y = [r], A = [Y], X = [], b = []

        let proof = create_mse_proof_ayxb(rng, &crs, vec![y], vec![r], vec![], vec![]);

        assert!(check_mse_proof_ayxb(&crs, vec![y], vec![], c, &proof));

        // Modify variable y
        let r_prime = Fr::rand(rng);
        let c_prime = c + y.mul(r_prime); // c' = c + r'Y = (r + r')Y

        let new_proof = proof.adapt_proof(&crs, vec![r_prime], vec![]);
        assert!(check_mse_proof_ayxb(
            &crs,
            vec![y],
            vec![],
            c_prime,
            &new_proof
        ));
    }
}
