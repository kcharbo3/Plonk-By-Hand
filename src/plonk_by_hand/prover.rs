use crate::field::Field;
use crate::math::ecc::{CurvePoint, ECC};
use crate::math::polynomial::Polynomial;
use crate::plonk_by_hand::constants;
use crate::plonk_by_hand::proof::{OpeningEvals, Proof, ProverPolys};
use crate::plonk_by_hand::public_coin::PublicCoin;
use crate::plonk_by_hand::structured_reference_string::SRS;
use crate::PythagoreanCircuit;

// TODO: make generic to any circuit, make circuit template/interface
pub struct Prover {
    py_circuit: PythagoreanCircuit,
    srs: SRS,
    ecc: ECC,
    field: Field,
    prover_polys: ProverPolys,
    pub opening_evals: OpeningEvals,
    pub_coin: PublicCoin,
}

impl Prover {
    pub fn new(field: Field, inputs: Vec<u32>, srs: SRS) -> Prover {
        let mut py_circuit = PythagoreanCircuit::new(field.clone());
        py_circuit.build_circuit_with_inputs(inputs);
        let ecc = ECC {
            field: Field { order: 101 },
        };

        Prover {
            py_circuit,
            srs,
            ecc,
            field,
            prover_polys: Default::default(),
            opening_evals: Default::default(),
            pub_coin: Default::default(),
        }
    }

    pub fn generate_proof(&mut self) -> Proof {
        self.set_wire_polys();
        self.set_z_poly();
        self.set_t_polys();

        self.set_first_opening_evals();

        self.set_r_poly();

        self.set_r_opening_eval();

        self.set_w_poly();
        self.set_wz_poly();

        Proof {
            a: self.commit_poly(&self.prover_polys.a),
            b: self.commit_poly(&self.prover_polys.b),
            c: self.commit_poly(&self.prover_polys.c),
            z: self.commit_poly(&self.prover_polys.z),
            t_lo: self.commit_poly(&self.prover_polys.t_lo),
            t_mid: self.commit_poly(&self.prover_polys.t_mid),
            t_hi: self.commit_poly(&self.prover_polys.t_hi),
            w: self.commit_poly(&self.prover_polys.w),
            wz: self.commit_poly(&self.prover_polys.wz),
            a_bar: self.opening_evals.a,
            b_bar: self.opening_evals.b,
            c_bar: self.opening_evals.c,
            left_copy_bar: self.opening_evals.left_copy,
            right_copy_bar: self.opening_evals.right_copy,
            r_bar: self.opening_evals.r,
            z_bar: self.opening_evals.z,
        }
    }

    pub fn set_public_coin(&mut self, public_coin: PublicCoin) {
        self.pub_coin = public_coin;
    }

    fn set_wire_polys(&mut self) {
        self.prover_polys.a = self.get_blinded_wire_poly(
            &self.py_circuit.circuit.circuit_polys.left,
            self.pub_coin.b2,
            self.pub_coin.b1,
        );

        self.prover_polys.b = self.get_blinded_wire_poly(
            &self.py_circuit.circuit.circuit_polys.right,
            self.pub_coin.b4,
            self.pub_coin.b3,
        );

        self.prover_polys.c = self.get_blinded_wire_poly(
            &self.py_circuit.circuit.circuit_polys.output,
            self.pub_coin.b6,
            self.pub_coin.b5,
        );
    }

    fn set_z_poly(&mut self) {
        self.prover_polys.z = self.get_blinded_z_poly(
            self.pub_coin.b9,
            self.pub_coin.b8,
            self.pub_coin.b7,
            self.pub_coin.beta,
            self.pub_coin.gamma,
        );
    }

    fn set_t_polys(&mut self) {
        let t = self.get_t_poly();
        let mut t_lo_coefficients = Vec::new();
        let mut t_mid_coefficients = Vec::new();
        let mut t_hi_coefficients = Vec::new();

        let num_coefficients = (t.degree + 1) / 3;
        for i in 0..num_coefficients {
            t_lo_coefficients.push(t.coefficients[i as usize]);
        }

        for i in 0..num_coefficients {
            t_mid_coefficients.push(t.coefficients[(num_coefficients + i) as usize]);
        }

        for i in 0..num_coefficients {
            t_hi_coefficients.push(t.coefficients[((num_coefficients * 2) + i) as usize]);
        }

        self.prover_polys.t = t.copy();
        self.prover_polys.t_lo = Polynomial {
            degree: num_coefficients - 1,
            coefficients: t_lo_coefficients,
            field: t.field.clone(),
        };

        self.prover_polys.t_mid = Polynomial {
            degree: num_coefficients - 1,
            coefficients: t_mid_coefficients,
            field: t.field.clone(),
        };

        self.prover_polys.t_hi = Polynomial {
            degree: num_coefficients - 1,
            coefficients: t_hi_coefficients,
            field: t.field,
        };
    }

    fn set_r_poly(&mut self) {
        self.prover_polys.r = self.get_r_poly();
    }

    fn set_w_poly(&mut self) {
        self.prover_polys.w = self.get_w_poly();
    }

    fn set_wz_poly(&mut self) {
        self.prover_polys.wz = self.get_wz_poly();
    }

    fn set_first_opening_evals(&mut self) {
        self.opening_evals.a = self.prover_polys.a.eval(self.pub_coin.zed);
        self.opening_evals.b = self.prover_polys.b.eval(self.pub_coin.zed);
        self.opening_evals.c = self.prover_polys.c.eval(self.pub_coin.zed);

        self.opening_evals.left_copy = self
            .py_circuit
            .circuit
            .circuit_polys
            .left_copy
            .eval(self.pub_coin.zed);
        self.opening_evals.right_copy = self
            .py_circuit
            .circuit
            .circuit_polys
            .right_copy
            .eval(self.pub_coin.zed);

        self.opening_evals.t = self.prover_polys.t.eval(self.pub_coin.zed);

        self.opening_evals.z = self.prover_polys.z.eval(
            self.field
                .multiply(self.pub_coin.zed, self.py_circuit.circuit.roots[1]),
        );
    }

    fn set_r_opening_eval(&mut self) {
        self.opening_evals.r = self.prover_polys.r.eval(self.pub_coin.zed);
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

    fn get_blinded_wire_poly(&self, wire_poly: &Polynomial, rand1: u32, rand2: u32) -> Polynomial {
        let blinded_z_poly = self.get_blinded_z_h_poly(1, vec![rand1, rand2]);

        Polynomial::poly_add(&blinded_z_poly, wire_poly)
    }

    fn get_blinded_z_poly(
        &mut self,
        rand1: u32,
        rand2: u32,
        rand3: u32,
        beta: u32,
        gamma: u32,
    ) -> Polynomial {
        let blinded_z_poly = self.get_blinded_z_h_poly(2, vec![rand1, rand2, rand3]);
        self.py_circuit.build_acc(beta, gamma);
        let acc_poly = &self.py_circuit.circuit.circuit_polys.acc;

        Polynomial::poly_add(&blinded_z_poly, acc_poly)
    }

    fn get_blinded_z_h_poly(&self, degree: u32, rands: Vec<u32>) -> Polynomial {
        let blinding_poly = Polynomial {
            degree,
            coefficients: rands,
            field: self.field.clone(),
        };

        Polynomial::poly_multiply(&self.py_circuit.circuit.circuit_polys.z_h, &blinding_poly)
    }

    fn get_t_poly(&self) -> Polynomial {
        let mut first_term = Polynomial::poly_multiply(&self.prover_polys.a, &self.prover_polys.b);
        first_term = Polynomial::poly_multiply(
            &first_term,
            &self.py_circuit.circuit.circuit_polys.multiply_selector,
        );
        first_term = Polynomial::poly_add(
            &first_term,
            &Polynomial::poly_multiply(
                &self.prover_polys.a,
                &self.py_circuit.circuit.circuit_polys.left_selector,
            ),
        );

        first_term = Polynomial::poly_add(
            &first_term,
            &Polynomial::poly_multiply(
                &self.prover_polys.b,
                &self.py_circuit.circuit.circuit_polys.right_selector,
            ),
        );

        first_term = Polynomial::poly_add(
            &first_term,
            &Polynomial::poly_multiply(
                &self.prover_polys.c,
                &self.py_circuit.circuit.circuit_polys.output_selector,
            ),
        );

        let helper1 = Polynomial {
            degree: 1,
            coefficients: vec![self.pub_coin.gamma, self.pub_coin.beta],
            field: self.field.clone(),
        };
        let helper2 = Polynomial {
            degree: 1,
            coefficients: vec![
                self.pub_coin.gamma,
                self.field.multiply(2, self.pub_coin.beta),
            ],
            field: self.field.clone(),
        };
        let helper3 = Polynomial {
            degree: 1,
            coefficients: vec![
                self.pub_coin.gamma,
                self.field.multiply(3, self.pub_coin.beta),
            ],
            field: self.field.clone(),
        };
        let mut second_term = Polynomial::poly_add(&self.prover_polys.a, &helper1);
        second_term = Polynomial::poly_multiply(
            &second_term,
            &Polynomial::poly_add(&self.prover_polys.b, &helper2),
        );
        second_term = Polynomial::poly_multiply(
            &second_term,
            &Polynomial::poly_add(&self.prover_polys.c, &helper3),
        );
        second_term = Polynomial::poly_multiply(&second_term, &self.prover_polys.z);
        second_term = Polynomial::scalar_multiply(&second_term, self.pub_coin.alpha);

        let gamma_poly = Polynomial {
            degree: 0,
            coefficients: vec![self.pub_coin.gamma],
            field: self.field.clone(),
        };
        let mut third_term = Polynomial::poly_add(
            &Polynomial::poly_add(
                &self.prover_polys.a,
                &self
                    .py_circuit
                    .circuit
                    .circuit_polys
                    .left_copy
                    .scalar_multiply(self.pub_coin.beta),
            ),
            &gamma_poly,
        );

        third_term = Polynomial::poly_multiply(
            &third_term,
            &Polynomial::poly_add(
                &self.prover_polys.b,
                &Polynomial::poly_add(
                    &self
                        .py_circuit
                        .circuit
                        .circuit_polys
                        .right_copy
                        .scalar_multiply(self.pub_coin.beta),
                    &gamma_poly,
                ),
            ),
        );

        third_term = Polynomial::poly_multiply(
            &third_term,
            &Polynomial::poly_add(
                &self.prover_polys.c,
                &Polynomial::poly_add(
                    &self
                        .py_circuit
                        .circuit
                        .circuit_polys
                        .output_copy
                        .scalar_multiply(self.pub_coin.beta),
                    &gamma_poly,
                ),
            ),
        );

        let mut z_at_w_x = Polynomial {
            degree: self.prover_polys.z.degree,
            coefficients: vec![0; self.prover_polys.z.coefficients.len()],
            field: self.prover_polys.z.field.clone(),
        };

        for i in 0..z_at_w_x.coefficients.len() {
            z_at_w_x.coefficients[i] = self.field.multiply(
                self.field
                    .exponent(self.py_circuit.circuit.roots[1], i as u32),
                self.prover_polys.z.coefficients[i],
            );
        }
        third_term = Polynomial::poly_multiply(&third_term, &z_at_w_x);
        third_term = third_term.scalar_multiply(self.pub_coin.alpha);
        third_term = Polynomial::additive_inverse_poly(&third_term);

        let inverse_one_poly = Polynomial::additive_inverse_poly(&Polynomial {
            degree: 0,
            coefficients: vec![1],
            field: self.field.clone(),
        });
        let mut fourth_term = Polynomial::poly_add(&self.prover_polys.z, &inverse_one_poly);
        fourth_term = Polynomial::poly_multiply(
            &fourth_term,
            &Polynomial::get_lagrange_1_poly(
                &self.field,
                (self.py_circuit.circuit.roots.len() - 1) as u32,
                &self.py_circuit.circuit.roots,
            ),
        );
        fourth_term = fourth_term.scalar_multiply(self.field.exponent(self.pub_coin.alpha, 2));

        let t_z = Polynomial::poly_add(
            &first_term,
            &Polynomial::poly_add(
                &second_term,
                &Polynomial::poly_add(&third_term, &fourth_term),
            ),
        );

        let (quotient, _remainder) =
            Polynomial::poly_divide(&t_z, &self.py_circuit.circuit.circuit_polys.z_h);

        quotient
    }

    fn get_r_poly(&self) -> Polynomial {
        let first_term = Polynomial::poly_add(
            &Polynomial::scalar_multiply(
                &self.py_circuit.circuit.circuit_polys.multiply_selector,
                self.field
                    .multiply(self.opening_evals.a, self.opening_evals.b),
            ),
            &Polynomial::poly_add(
                &Polynomial::scalar_multiply(
                    &self.py_circuit.circuit.circuit_polys.left_selector,
                    self.opening_evals.a,
                ),
                &Polynomial::poly_add(
                    &Polynomial::scalar_multiply(
                        &self.py_circuit.circuit.circuit_polys.right_selector,
                        self.opening_evals.b,
                    ),
                    &Polynomial::scalar_multiply(
                        &self.py_circuit.circuit.circuit_polys.output_selector,
                        self.opening_evals.c,
                    ),
                ),
            ),
        );

        let second_term_scalar = self.field.multiply(
            self.field.multiply(
                self.field.add(
                    self.opening_evals.a,
                    self.field.add(
                        self.field.multiply(self.pub_coin.beta, self.pub_coin.zed),
                        self.pub_coin.gamma,
                    ),
                ),
                self.field.multiply(
                    self.field.add(
                        self.opening_evals.b,
                        self.field.add(
                            self.field.multiply(
                                self.field.multiply(self.pub_coin.beta, 2),
                                self.pub_coin.zed,
                            ),
                            self.pub_coin.gamma,
                        ),
                    ),
                    self.field.add(
                        self.opening_evals.c,
                        self.field.add(
                            self.field.multiply(
                                self.field.multiply(self.pub_coin.beta, 3),
                                self.pub_coin.zed,
                            ),
                            self.pub_coin.gamma,
                        ),
                    ),
                ),
            ),
            self.pub_coin.alpha,
        );

        let second_term = Polynomial::scalar_multiply(&self.prover_polys.z, second_term_scalar);

        let third_term_scalar = self.field.multiply(
            self.field.multiply(
                self.field.add(
                    self.opening_evals.a,
                    self.field.add(
                        self.field
                            .multiply(self.pub_coin.beta, self.opening_evals.left_copy),
                        self.pub_coin.gamma,
                    ),
                ),
                self.field.multiply(
                    self.field.add(
                        self.opening_evals.b,
                        self.field.add(
                            self.field
                                .multiply(self.pub_coin.beta, self.opening_evals.right_copy),
                            self.pub_coin.gamma,
                        ),
                    ),
                    self.field
                        .multiply(self.pub_coin.beta, self.opening_evals.z),
                ),
            ),
            self.pub_coin.alpha,
        );
        let third_term = Polynomial::additive_inverse_poly(&Polynomial::scalar_multiply(
            &self.py_circuit.circuit.circuit_polys.output_copy,
            third_term_scalar,
        ));

        let fourth_term_scalar = self.field.multiply(
            Polynomial::get_lagrange_1_poly(
                &self.field,
                (self.py_circuit.circuit.roots.len() - 1) as u32,
                &self.py_circuit.circuit.roots,
            )
            .eval(self.pub_coin.zed),
            self.field.exponent(self.pub_coin.alpha, 2),
        );
        let fourth_term = Polynomial::scalar_multiply(&self.prover_polys.z, fourth_term_scalar);

        Polynomial::poly_add(
            &first_term,
            &Polynomial::poly_add(
                &second_term,
                &Polynomial::poly_add(&third_term, &fourth_term),
            ),
        )
    }

    fn get_w_poly(&self) -> Polynomial {
        let mut first_term = Polynomial::poly_add(
            &self.prover_polys.t_lo,
            &Polynomial::poly_add(
                &Polynomial::scalar_multiply(
                    &self.prover_polys.t_mid,
                    self.field.exponent(
                        self.pub_coin.zed,
                        (self.py_circuit.circuit.roots.len() + 2) as u32,
                    ),
                ),
                &Polynomial::scalar_multiply(
                    &self.prover_polys.t_hi,
                    self.field.exponent(
                        self.pub_coin.zed,
                        ((2 * self.py_circuit.circuit.roots.len()) + 4) as u32,
                    ),
                ),
            ),
        );

        first_term = Polynomial::poly_add(
            &first_term,
            &Polynomial {
                degree: 0,
                coefficients: vec![self.field.additive_inverse(self.opening_evals.t)],
                field: self.field.clone(),
            },
        );

        let mut second_term = Polynomial::scalar_multiply(
            &Polynomial::poly_add(
                &self.prover_polys.r,
                &Polynomial::from_scalar(
                    self.field.additive_inverse(self.opening_evals.r),
                    self.field.clone(),
                ),
            ),
            self.pub_coin.v,
        );

        second_term = Polynomial::poly_add(
            &second_term,
            &Polynomial::scalar_multiply(
                &Polynomial::poly_add(
                    &self.prover_polys.a,
                    &Polynomial::from_scalar(
                        self.field.additive_inverse(self.opening_evals.a),
                        self.field.clone(),
                    ),
                ),
                self.field.exponent(self.pub_coin.v, 2),
            ),
        );

        second_term = Polynomial::poly_add(
            &second_term,
            &Polynomial::scalar_multiply(
                &Polynomial::poly_add(
                    &self.prover_polys.b,
                    &Polynomial::from_scalar(
                        self.field.additive_inverse(self.opening_evals.b),
                        self.field.clone(),
                    ),
                ),
                self.field.exponent(self.pub_coin.v, 3),
            ),
        );

        second_term = Polynomial::poly_add(
            &second_term,
            &Polynomial::scalar_multiply(
                &Polynomial::poly_add(
                    &self.prover_polys.c,
                    &Polynomial::from_scalar(
                        self.field.additive_inverse(self.opening_evals.c),
                        self.field.clone(),
                    ),
                ),
                self.field.exponent(self.pub_coin.v, 4),
            ),
        );

        second_term = Polynomial::poly_add(
            &second_term,
            &Polynomial::scalar_multiply(
                &Polynomial::poly_add(
                    &self.py_circuit.circuit.circuit_polys.left_copy,
                    &Polynomial::from_scalar(
                        self.field.additive_inverse(self.opening_evals.left_copy),
                        self.field.clone(),
                    ),
                ),
                self.field.exponent(self.pub_coin.v, 5),
            ),
        );

        second_term = Polynomial::poly_add(
            &second_term,
            &Polynomial::scalar_multiply(
                &Polynomial::poly_add(
                    &self.py_circuit.circuit.circuit_polys.right_copy,
                    &Polynomial::from_scalar(
                        self.field.additive_inverse(self.opening_evals.right_copy),
                        self.field.clone(),
                    ),
                ),
                self.field.exponent(self.pub_coin.v, 6),
            ),
        );

        let denominator = Polynomial {
            degree: 1,
            coefficients: vec![self.field.additive_inverse(self.pub_coin.zed), 1],
            field: self.field.clone(),
        };

        let (quotient, _remainder) = Polynomial::poly_divide(
            &Polynomial::poly_add(&first_term, &second_term),
            &denominator,
        );
        quotient
    }

    fn get_wz_poly(&self) -> Polynomial {
        let numerator = Polynomial::poly_add(
            &self.prover_polys.z,
            &Polynomial::from_scalar(
                self.field.additive_inverse(self.opening_evals.z),
                self.field.clone(),
            ),
        );
        let denominator = Polynomial {
            degree: 1,
            coefficients: vec![
                self.field.additive_inverse(
                    self.field
                        .multiply(self.pub_coin.zed, self.py_circuit.circuit.roots[1]),
                ),
                1,
            ],
            field: self.field.clone(),
        };

        let (quotient, _remainder) = Polynomial::poly_divide(&numerator, &denominator);
        quotient
    }
}

fn test_setup_prover() -> Prover {
    let mut srs = constants::SRS_BY_HAND;
    srs.generate_g_1_points();

    let mut prover = Prover::new(constants::FIELD_17.clone(), vec![3, 4, 5], srs);
    prover.set_public_coin(constants::PUB_COIN.clone());

    prover.generate_proof();

    prover
}

#[test]
fn test_blinded_wire_polys() {
    let prover = test_setup_prover();

    let a_coefficients = prover
        .get_blinded_wire_poly(&prover.py_circuit.circuit.circuit_polys.left, 4, 7)
        .coefficients;
    assert_eq!(a_coefficients, vec![14, 6, 3, 3, 4, 7]);

    let b_coefficients = prover
        .get_blinded_wire_poly(&prover.py_circuit.circuit.circuit_polys.right, 12, 11)
        .coefficients;
    assert_eq!(b_coefficients, vec![12, 9, 14, 13, 12, 11]);

    let c_coefficients = prover
        .get_blinded_wire_poly(&prover.py_circuit.circuit.circuit_polys.output, 2, 16)
        .coefficients;
    assert_eq!(c_coefficients, vec![4, 6, 11, 4, 2, 16]);

    let a_commitment = prover.commit_poly(&prover.prover_polys.a);
    assert!(CurvePoint::equals(
        &a_commitment,
        &CurvePoint { x: 91, y: 66 }
    ));

    let b_commitment = prover.commit_poly(&prover.prover_polys.b);
    assert!(CurvePoint::equals(
        &b_commitment,
        &CurvePoint { x: 26, y: 45 }
    ));

    let c_commitment = prover.commit_poly(&prover.prover_polys.c);
    assert!(CurvePoint::equals(
        &c_commitment,
        &CurvePoint { x: 91, y: 35 }
    ));
}

#[test]
fn test_accumulator() {
    let prover = test_setup_prover();

    assert!(CurvePoint::equals(
        &prover.commit_poly(&prover.prover_polys.z),
        &CurvePoint { x: 32, y: 59 }
    ));
}

#[test]
fn test_t_poly() {
    let prover = test_setup_prover();

    let t_poly = prover.get_t_poly();
    assert_eq!(
        t_poly.coefficients,
        vec![11, 16, 13, 9, 0, 13, 13, 8, 1, 2, 10, 1, 15, 6, 16, 2, 7, 11]
    );

    assert_eq!(
        prover.prover_polys.t_lo.coefficients,
        vec![11, 16, 13, 9, 0, 13]
    );
    assert_eq!(
        prover.prover_polys.t_mid.coefficients,
        vec![13, 8, 1, 2, 10, 1]
    );
    assert_eq!(
        prover.prover_polys.t_hi.coefficients,
        vec![15, 6, 16, 2, 7, 11]
    );

    let t_lo_com = prover.commit_poly(&prover.prover_polys.t_lo);
    let t_mid_com = prover.commit_poly(&prover.prover_polys.t_mid);
    let t_hi_com = prover.commit_poly(&prover.prover_polys.t_hi);

    assert!(CurvePoint::equals(&t_lo_com, &CurvePoint { x: 12, y: 32 }));
    assert!(CurvePoint::equals(&t_mid_com, &CurvePoint { x: 26, y: 45 }));
    assert!(CurvePoint::equals(&t_hi_com, &CurvePoint { x: 91, y: 66 }));
}

#[test]
fn test_opening_evals() {
    let prover = test_setup_prover();

    let actual_evals = vec![
        prover.opening_evals.a,
        prover.opening_evals.b,
        prover.opening_evals.c,
        prover.opening_evals.left_copy,
        prover.opening_evals.right_copy,
        prover.opening_evals.t,
        prover.opening_evals.z,
    ];
    let expected_evals = vec![15, 13, 5, 1, 12, 1, 15];
    assert_eq!(actual_evals, expected_evals);

    let r_eval = prover.opening_evals.r;
    let expected_r_eval = 15;
    assert_eq!(r_eval, expected_r_eval);
}

#[test]
fn test_r_poly() {
    let prover = test_setup_prover();

    let actual_coefficients = prover.prover_polys.r.coefficients;
    let expected_coefficients = vec![0, 16, 9, 13, 8, 15, 16];
    assert_eq!(actual_coefficients, expected_coefficients);
}

#[test]
fn test_w_polys() {
    let prover = test_setup_prover();

    let actual_coefficients = prover.prover_polys.w.coefficients.clone();
    let expected_coefficients = vec![16, 13, 2, 9, 3, 5];
    assert_eq!(actual_coefficients, expected_coefficients);

    let actual_coefficients = prover.prover_polys.wz.coefficients;
    let expected_coefficients = vec![13, 14, 2, 13, 2, 14];
    assert_eq!(actual_coefficients, expected_coefficients);
}

#[test]
fn test_w_commits() {
    let prover = test_setup_prover();

    let w_com = prover.commit_poly(&prover.prover_polys.w);
    let wz_com = prover.commit_poly(&prover.prover_polys.wz);

    assert!(CurvePoint::equals(&w_com, &CurvePoint { x: 91, y: 35 }));
    assert!(CurvePoint::equals(&wz_com, &CurvePoint { x: 65, y: 98 }));
}
