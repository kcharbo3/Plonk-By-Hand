use crate::field::Field;
use crate::math::ecc::{ComplexScalar, ECC};
use crate::math::pairing;
use crate::math::polynomial::Polynomial;
use crate::plonk_by_hand::constants;
use crate::plonk_by_hand::proof::Proof;
use crate::plonk_by_hand::public_coin::PublicCoin;
use crate::plonk_by_hand::structured_reference_string::SRS;
use crate::{CurvePoint, Prover, PythagoreanCircuit};

pub struct Verifier {
    pub circuit: PythagoreanCircuit,
    pub srs: SRS,
    pub ecc: ECC,
    pub field: Field,
    pub pub_coin: PublicCoin,
    pub proof: Proof,
    pub commitments: VerifierCommitments,
    pub vals: VerifierVals,
}

#[derive(Default)]
pub struct VerifierCommitments {
    pub left_selector: CurvePoint,
    pub right_selector: CurvePoint,
    pub output_selector: CurvePoint,
    pub multiply_selector: CurvePoint,
    pub c_selector: CurvePoint,
    pub left_copy: CurvePoint,
    pub right_copy: CurvePoint,
    pub output_copy: CurvePoint,
}

#[derive(Default)]
pub struct VerifierVals {
    pub z_h_opening: u32,
    pub lagrange_1_opening: u32,
    pub t_opening: u32,
    pub d_commitment: CurvePoint,
    pub f_commitment: CurvePoint,
    pub e_commitment: CurvePoint,
}

impl Verifier {
    pub fn new(field: Field, srs: SRS) -> Verifier {
        let mut py_circuit = PythagoreanCircuit::new(field.clone());
        py_circuit.build_circuit();
        let ecc = ECC {
            field: Field { order: 101 },
        };

        Verifier {
            circuit: py_circuit,
            srs,
            ecc,
            field,
            pub_coin: Default::default(),
            proof: Default::default(),
            commitments: Default::default(),
            vals: Default::default(),
        }
    }

    pub fn verify_proof(&mut self) -> bool {
        // Step 1
        let mut verified = self.verify_commitments_on_curve();

        // Step 2
        verified = verified && self.verify_openings_in_field();

        // Skip step 3 since no public inputs

        // Step 4
        self.set_z_h_opening();

        // Step 5
        self.set_lagrange_1_opening();

        // Skip step 6 since no public inputs

        // Step 7
        self.set_t_opening();

        // Step 8
        self.set_d_commitment();

        // Step 9
        self.set_f_commitment();

        // Step 10
        self.set_e_commitment();

        // Step 11
        verified && self.check_pairing()
    }

    pub fn provide_proof(&mut self, pub_coin: PublicCoin, proof: Proof) {
        self.pub_coin = pub_coin;
        self.proof = proof;
    }

    pub fn preprocess(&mut self) {
        self.commitments = VerifierCommitments {
            left_selector: self.commit_poly(&self.circuit.circuit.circuit_polys.left_selector),
            right_selector: self.commit_poly(&self.circuit.circuit.circuit_polys.right_selector),
            output_selector: self.commit_poly(&self.circuit.circuit.circuit_polys.output_selector),
            multiply_selector: self
                .commit_poly(&self.circuit.circuit.circuit_polys.multiply_selector),
            c_selector: CurvePoint::point_at_infinity(),
            left_copy: self.commit_poly(&self.circuit.circuit.circuit_polys.left_copy),
            right_copy: self.commit_poly(&self.circuit.circuit.circuit_polys.right_copy),
            output_copy: self.commit_poly(&self.circuit.circuit.circuit_polys.output_copy),
        }
    }

    pub fn verify_commitments_on_curve(&self) -> bool {
        self.in_curve(&self.proof.a)
            && self.in_curve(&self.proof.b)
            && self.in_curve(&self.proof.c)
            && self.in_curve(&self.proof.z)
            && self.in_curve(&self.proof.t_lo)
            && self.in_curve(&self.proof.t_mid)
            && self.in_curve(&self.proof.t_hi)
            && self.in_curve(&self.proof.w)
            && self.in_curve(&self.proof.wz)
    }

    pub fn verify_openings_in_field(&self) -> bool {
        self.in_scalar_field(self.proof.a_bar)
            && self.in_scalar_field(self.proof.b_bar)
            && self.in_scalar_field(self.proof.c_bar)
            && self.in_scalar_field(self.proof.left_copy_bar)
            && self.in_scalar_field(self.proof.right_copy_bar)
            && self.in_scalar_field(self.proof.r_bar)
            && self.in_scalar_field(self.proof.z_bar)
    }

    pub fn set_z_h_opening(&mut self) {
        self.vals.z_h_opening = self
            .circuit
            .circuit
            .circuit_polys
            .z_h
            .eval(self.pub_coin.zed);
    }

    pub fn set_lagrange_1_opening(&mut self) {
        self.vals.lagrange_1_opening = Polynomial::get_lagrange_1_poly(
            &self.field,
            (self.circuit.circuit.roots.len() - 1) as u32,
            &self.circuit.circuit.roots,
        )
        .eval(self.pub_coin.zed);
    }

    pub fn set_t_opening(&mut self) {
        let first_term = self.field.add(
            self.proof.a_bar,
            self.field.add(
                self.field
                    .multiply(self.pub_coin.beta, self.proof.left_copy_bar),
                self.pub_coin.gamma,
            ),
        );

        let second_term = self.field.add(
            self.proof.b_bar,
            self.field.add(
                self.field
                    .multiply(self.pub_coin.beta, self.proof.right_copy_bar),
                self.pub_coin.gamma,
            ),
        );

        let third_term = self.field.add(self.proof.c_bar, self.pub_coin.gamma);

        let middle_term = self.field.multiply(
            first_term,
            self.field.multiply(
                second_term,
                self.field.multiply(
                    third_term,
                    self.field.multiply(self.proof.z_bar, self.pub_coin.alpha),
                ),
            ),
        );

        let mut res = self
            .field
            .add(self.proof.r_bar, self.field.additive_inverse(middle_term));
        res = self.field.add(
            res,
            self.field.additive_inverse(self.field.multiply(
                self.vals.lagrange_1_opening,
                self.field.exponent(self.pub_coin.alpha, 2),
            )),
        );

        res = self.field.divide(res, self.vals.z_h_opening);

        self.vals.t_opening = res;
    }

    pub fn set_d_commitment(&mut self) {
        let first_term = self.ecc.multiply(
            self.field.multiply(
                self.proof.a_bar,
                self.field.multiply(self.proof.b_bar, self.pub_coin.v),
            ),
            &self.commitments.multiply_selector,
        );

        let second_term = self.ecc.multiply(
            self.field.multiply(self.proof.a_bar, self.pub_coin.v),
            &self.commitments.left_selector,
        );

        let third_term = self.ecc.multiply(
            self.field.multiply(self.proof.b_bar, self.pub_coin.v),
            &self.commitments.right_selector,
        );

        let fourth_term = self.ecc.multiply(
            self.field.multiply(self.proof.c_bar, self.pub_coin.v),
            &self.commitments.output_selector,
        );

        let fifth_term = self
            .ecc
            .multiply(self.pub_coin.v, &self.commitments.c_selector);

        let mut res = self.ecc.add(
            &first_term,
            &self.ecc.add(
                &second_term,
                &self
                    .ecc
                    .add(&third_term, &self.ecc.add(&fourth_term, &fifth_term)),
            ),
        );

        let mut middle_term_scalar = self.field.add(
            self.proof.a_bar,
            self.field.add(
                self.field.multiply(self.pub_coin.beta, self.pub_coin.zed),
                self.pub_coin.gamma,
            ),
        );
        middle_term_scalar = self.field.multiply(
            middle_term_scalar,
            self.field.add(
                self.proof.b_bar,
                self.field.add(
                    self.field.multiply(
                        self.field.multiply(self.pub_coin.beta, 2),
                        self.pub_coin.zed,
                    ),
                    self.pub_coin.gamma,
                ),
            ),
        );
        middle_term_scalar = self.field.multiply(
            middle_term_scalar,
            self.field.add(
                self.proof.c_bar,
                self.field.add(
                    self.field.multiply(
                        self.field.multiply(self.pub_coin.beta, 3),
                        self.pub_coin.zed,
                    ),
                    self.pub_coin.gamma,
                ),
            ),
        );
        middle_term_scalar = self.field.multiply(
            middle_term_scalar,
            self.field.multiply(self.pub_coin.alpha, self.pub_coin.v),
        );
        middle_term_scalar = self.field.add(
            middle_term_scalar,
            self.field.multiply(
                self.vals.lagrange_1_opening,
                self.field
                    .multiply(self.field.exponent(self.pub_coin.alpha, 2), self.pub_coin.v),
            ),
        );

        middle_term_scalar = self.field.add(middle_term_scalar, self.pub_coin.u);
        let middle_term = self.ecc.multiply(middle_term_scalar, &self.proof.z);

        res = self.ecc.add(&res, &middle_term);

        let mut last_term_scalar = self.field.add(
            self.proof.a_bar,
            self.field.add(
                self.field
                    .multiply(self.pub_coin.beta, self.proof.left_copy_bar),
                self.pub_coin.gamma,
            ),
        );
        last_term_scalar = self.field.multiply(
            last_term_scalar,
            self.field.add(
                self.proof.b_bar,
                self.field.add(
                    self.field
                        .multiply(self.pub_coin.beta, self.proof.right_copy_bar),
                    self.pub_coin.gamma,
                ),
            ),
        );
        last_term_scalar = self.field.multiply(
            last_term_scalar,
            self.field.multiply(
                self.pub_coin.alpha,
                self.field.multiply(
                    self.pub_coin.v,
                    self.field.multiply(self.pub_coin.beta, self.proof.z_bar),
                ),
            ),
        );
        last_term_scalar = self.field.additive_inverse(last_term_scalar);
        let last_term = self
            .ecc
            .multiply(last_term_scalar, &self.commitments.output_copy);

        res = self.ecc.add(&res, &last_term);
        self.vals.d_commitment = res;
    }

    pub fn set_f_commitment(&mut self) {
        let mut first_term = self.ecc.add(
            &self.proof.t_lo,
            &self.ecc.multiply(
                self.field.exponent(
                    self.pub_coin.zed,
                    (self.circuit.circuit.roots.len() + 2) as u32,
                ),
                &self.proof.t_mid,
            ),
        );
        first_term = self.ecc.add(
            &first_term,
            &self.ecc.multiply(
                self.field.exponent(
                    self.pub_coin.zed,
                    self.field
                        .multiply(self.circuit.circuit.roots.len() as u32, 2)
                        + 4,
                ),
                &self.proof.t_hi,
            ),
        );

        let second_term = self.vals.d_commitment.clone();

        let third_term = self.ecc.add(
            &self
                .ecc
                .multiply(self.field.exponent(self.pub_coin.v, 2), &self.proof.a),
            &self.ecc.add(
                &self
                    .ecc
                    .multiply(self.field.exponent(self.pub_coin.v, 3), &self.proof.b),
                &self.ecc.add(
                    &self
                        .ecc
                        .multiply(self.field.exponent(self.pub_coin.v, 4), &self.proof.c),
                    &self.ecc.add(
                        &self.ecc.multiply(
                            self.field.exponent(self.pub_coin.v, 5),
                            &self.commitments.left_copy,
                        ),
                        &self.ecc.multiply(
                            self.field.exponent(self.pub_coin.v, 6),
                            &self.commitments.right_copy,
                        ),
                    ),
                ),
            ),
        );

        self.vals.f_commitment = self
            .ecc
            .add(&first_term, &self.ecc.add(&second_term, &third_term));
    }

    pub fn set_e_commitment(&mut self) {
        let scalar_term = self.field.add(
            self.vals.t_opening,
            self.field.add(
                self.field.multiply(self.pub_coin.v, self.proof.r_bar),
                self.field.add(
                    self.field
                        .multiply(self.field.exponent(self.pub_coin.v, 2), self.proof.a_bar),
                    self.field.add(
                        self.field
                            .multiply(self.field.exponent(self.pub_coin.v, 3), self.proof.b_bar),
                        self.field.add(
                            self.field.multiply(
                                self.field.exponent(self.pub_coin.v, 4),
                                self.proof.c_bar,
                            ),
                            self.field.add(
                                self.field.multiply(
                                    self.field.exponent(self.pub_coin.v, 5),
                                    self.proof.left_copy_bar,
                                ),
                                self.field.add(
                                    self.field.multiply(
                                        self.field.exponent(self.pub_coin.v, 6),
                                        self.proof.right_copy_bar,
                                    ),
                                    self.field.multiply(self.pub_coin.u, self.proof.z_bar),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        );

        self.vals.e_commitment = self.ecc.multiply(scalar_term, &self.srs.g_1);
    }

    pub fn check_pairing(&mut self) -> bool {
        let left_first_term = self.ecc.add(
            &self.proof.w,
            &self.ecc.multiply(self.pub_coin.u, &self.proof.wz),
        );

        let mut right_first_term = self.ecc.add(
            &self.ecc.multiply(self.pub_coin.zed, &self.proof.w),
            &self.ecc.multiply(
                self.field.multiply(
                    self.pub_coin.u,
                    self.field
                        .multiply(self.pub_coin.zed, self.circuit.circuit.roots[1]),
                ),
                &self.proof.wz,
            ),
        );

        right_first_term = self.ecc.add(
            &right_first_term,
            &self.ecc.add(
                &self.vals.f_commitment,
                &self.ecc.inversion(&self.vals.e_commitment),
            ),
        );

        let left_multiple = self.ecc.get_generator_multiple(
            &left_first_term,
            &self.srs.g_1_points[0],
            self.field.order,
        );
        let right_multiple = self.ecc.get_generator_multiple(
            &right_first_term,
            &self.srs.g_1_points[0],
            self.field.order,
        );

        let pairing = pairing::Pairing {
            ecc: self.ecc.clone(),
            r: self.field.order,
        };

        let mut left_pairing =
            pairing.get_base_pairing(&self.srs.g_2_points[1], &self.srs.g_1_points[0]);
        let mut right_pairing =
            pairing.get_base_pairing(&self.srs.g_2_points[0], &self.srs.g_1_points[0]);

        left_pairing = ComplexScalar::exponent(&self.ecc.field, &left_pairing, left_multiple);
        right_pairing = ComplexScalar::exponent(&self.ecc.field, &right_pairing, right_multiple);

        ComplexScalar::equals(&left_pairing, &right_pairing)
    }

    fn in_scalar_field(&self, scalar: u32) -> bool {
        scalar < self.field.order
    }

    fn in_curve(&self, point: &CurvePoint) -> bool {
        if point.x >= self.ecc.field.order || point.y >= self.ecc.field.order {
            return false;
        }
        // must satisfy y^2 = x^3 + 3
        let left = self.ecc.field.exponent(point.y, 2);
        let right = self.ecc.field.add(self.ecc.field.exponent(point.x, 3), 3);

        left == right
    }

    fn commit_poly(&self, poly: &Polynomial) -> CurvePoint {
        let mut res = CurvePoint::point_at_infinity();

        for i in 0..poly.coefficients.len() {
            res = self.ecc.add(
                &res,
                &self
                    .ecc
                    .multiply(poly.coefficients[i], &self.srs.g_1_points[i]),
            );
        }

        res
    }
}

fn test_setup_verifier_with_proof() -> Verifier {
    let field_17 = constants::FIELD_17.clone();

    let mut srs = constants::SRS_BY_HAND.clone();
    srs.generate_g_1_points();
    srs.generate_g_2_points();

    let pub_coin = constants::PUB_COIN.clone();

    const A: u32 = 3;
    const B: u32 = 4;
    const C: u32 = 5;

    let mut prover = Prover::new(field_17.clone(), vec![A, B, C], srs.copy());
    prover.set_public_coin(pub_coin.clone());

    let mut verifier = Verifier::new(field_17.clone(), srs.copy());
    verifier.preprocess();
    verifier.provide_proof(pub_coin.clone(), prover.generate_proof());
    verifier
}

#[test]
fn test_preprocess() {
    let verifier = test_setup_verifier_with_proof();

    assert!(CurvePoint::equals(
        &verifier.commitments.left_selector,
        &CurvePoint::new(32, 42)
    ));
    assert!(CurvePoint::equals(
        &verifier.commitments.right_selector,
        &CurvePoint::new(32, 42)
    ));
    assert!(CurvePoint::equals(
        &verifier.commitments.output_selector,
        &CurvePoint::new(1, 99)
    ));
    assert!(CurvePoint::equals(
        &verifier.commitments.multiply_selector,
        &CurvePoint::new(12, 69)
    ));
    assert!(CurvePoint::equals(
        &verifier.commitments.c_selector,
        &CurvePoint::point_at_infinity()
    ));
    assert!(CurvePoint::equals(
        &verifier.commitments.left_copy,
        &CurvePoint::new(68, 74)
    ));
    assert!(CurvePoint::equals(
        &verifier.commitments.right_copy,
        &CurvePoint::new(65, 3)
    ));
    assert!(CurvePoint::equals(
        &verifier.commitments.output_copy,
        &CurvePoint::new(18, 49)
    ));
}

#[test]
fn test_verify_commitments_in_field() {
    let verifier = test_setup_verifier_with_proof();

    assert!(verifier.verify_commitments_on_curve());
}

#[test]
fn test_verify_openings_in_field() {
    let verifier = test_setup_verifier_with_proof();

    assert!(verifier.verify_openings_in_field());
}

#[test]
fn test_set_z_h_and_lagrange_opening() {
    let mut verifier = test_setup_verifier_with_proof();
    verifier.set_z_h_opening();
    verifier.set_lagrange_1_opening();
    assert_eq!(12, verifier.vals.z_h_opening);
    assert_eq!(5, verifier.vals.lagrange_1_opening);
}

#[test]
fn test_set_t_opening() {
    let mut verifier = test_setup_verifier_with_proof();
    verifier.set_z_h_opening();
    verifier.set_lagrange_1_opening();
    verifier.set_t_opening();
    assert_eq!(1, verifier.vals.t_opening);
}

#[test]
fn test_d_commitment() {
    let mut verifier = test_setup_verifier_with_proof();
    verifier.set_z_h_opening();
    verifier.set_lagrange_1_opening();
    verifier.set_t_opening();
    verifier.set_d_commitment();
    assert!(CurvePoint::equals(
        &CurvePoint::point_at_infinity(),
        &verifier.vals.d_commitment
    ));
}

#[test]
fn test_f_commitment() {
    let mut verifier = test_setup_verifier_with_proof();
    verifier.set_z_h_opening();
    verifier.set_lagrange_1_opening();
    verifier.set_t_opening();
    verifier.set_d_commitment();
    verifier.set_f_commitment();
    assert!(CurvePoint::equals(
        &CurvePoint::new(68, 27),
        &verifier.vals.f_commitment
    ));
}

#[test]
fn test_e_commitment() {
    let mut verifier = test_setup_verifier_with_proof();
    verifier.set_z_h_opening();
    verifier.set_lagrange_1_opening();
    verifier.set_t_opening();
    verifier.set_d_commitment();
    verifier.set_f_commitment();
    verifier.set_e_commitment();
    assert!(CurvePoint::equals(
        &CurvePoint::new(1, 2),
        &verifier.vals.e_commitment
    ));
}

#[test]
fn test_pairing() {
    let mut verifier = test_setup_verifier_with_proof();
    verifier.set_z_h_opening();
    verifier.set_lagrange_1_opening();
    verifier.set_t_opening();
    verifier.set_d_commitment();
    verifier.set_f_commitment();
    verifier.set_e_commitment();

    assert!(verifier.check_pairing());
}
