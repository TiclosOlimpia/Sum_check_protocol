use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial};

mod verifier;
mod prover;

use verifier::Fq;

fn main() {
    let v = 3;
    // Create a multivariate polynomial in 3 variables:
    // f(x_0, x_1, x_2) = 2*x_0^3 + x_0*x_2 + x_1*x_2 
    let g: prover::MultiPoly = SparsePolynomial::from_coefficients_vec(
        v,
        vec![
            (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
            (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
            (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),  
        ],
    );

    let r_random_vector = vec![2,3,6];

    let mut v = verifier::Verifier::new(&g, r_random_vector);
    v.verify_sum_check_protocol();
}