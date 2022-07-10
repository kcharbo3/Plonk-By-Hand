use crate::math::field::Field;
use crate::math::polynomial::{Point, Polynomial};

// user must: set_inputs, insert_gate (s), compute_witness, then build_polynomials
pub struct Circuit {
    pub gates: Vec<Gate>,
    pub inputs: Vec<(u32, usize)>, // (input val, index in witness)
    pub witness: Vec<u32>,
    pub field: Field,
    pub coset1: Vec<u32>,
    pub coset2: Vec<u32>,
    pub roots: Vec<u32>,
    pub circuit_polys: CircuitPolynomials,
}

impl Circuit {
    pub fn insert_gate(&mut self, gate: Gate) {
        self.gates.push(gate);
    }

    pub fn set_inputs(&mut self, inputs: Vec<(u32, usize)>) {
        self.inputs = inputs;
        for (val, index) in self.inputs.iter() {
            self.witness[*index] = *val;
        }
    }

    pub fn compute_witness(&mut self) {
        for gate in self.gates.iter() {
            let left_val = self.witness[gate.left_index];
            let right_val = self.witness[gate.right_index];

            match gate {
                Gate {
                    gate_type: GateType::Addition,
                    ..
                } => {
                    self.witness[gate.output_index] = self.field.add(left_val, right_val);
                }
                Gate {
                    gate_type: GateType::Multiplication,
                    ..
                } => {
                    self.witness[gate.output_index] = self.field.multiply(left_val, right_val);
                }
            }
        }
    }

    pub fn build_polynomials_with_input(&mut self) {
        self.circuit_polys.left = self.build_poly(self.get_left_inputs());
        self.circuit_polys.right = self.build_poly(self.get_right_inputs());
        self.circuit_polys.output = self.build_poly(self.get_outputs());
    }

    pub fn build_polynomials(&mut self) {
        self.circuit_polys.left_selector = self.build_poly(self.get_left_selectors());
        self.circuit_polys.right_selector = self.build_poly(self.get_right_selectors());
        self.circuit_polys.output_selector = self.build_poly(self.get_output_selectors());
        self.circuit_polys.multiply_selector = self.build_poly(self.get_multiply_selectors());
        self.circuit_polys.left_copy = self.build_poly(self.get_left_copy_constraints());
        self.circuit_polys.right_copy = self.build_poly(self.get_right_copy_constraints());
        self.circuit_polys.output_copy = self.build_poly(self.get_output_copy_constraints());
        self.circuit_polys.z_h = self.get_z_h();
    }

    pub fn build_accumulator(&mut self, beta: u32, gamma: u32) {
        // left, right, and output copy polys must be initialized first
        self.circuit_polys.acc = self.build_poly(self.get_accumulator_values(
            self.circuit_polys.left.degree,
            beta,
            gamma,
        ));
    }

    pub fn build_poly(&self, y_vals: Vec<u32>) -> Polynomial {
        let points = y_vals
            .iter()
            .zip(self.roots.iter())
            .map(|(input, root)| Point {
                x: *root,
                y: *input,
            })
            .collect();

        Polynomial::create_poly_from_points(points, &self.field)
    }

    pub fn get_z_h(&self) -> Polynomial {
        let degree = self.roots.len() as u32;
        let mut coefficients = vec![0; (degree + 1) as usize];
        coefficients[0] = self.field.additive_inverse(1);
        coefficients[degree as usize] = 1;
        Polynomial {
            degree,
            coefficients,
            field: self.field.clone(),
        }
    }

    pub fn get_left_inputs(&self) -> Vec<u32> {
        let mut left_inputs: Vec<u32> = Vec::new();
        for gate in self.gates.iter() {
            left_inputs.push(self.witness[gate.left_index]);
        }

        left_inputs
    }

    pub fn get_right_inputs(&self) -> Vec<u32> {
        let mut right_inputs: Vec<u32> = Vec::new();
        for gate in self.gates.iter() {
            right_inputs.push(self.witness[gate.right_index]);
        }

        right_inputs
    }

    pub fn get_outputs(&self) -> Vec<u32> {
        let mut outputs: Vec<u32> = Vec::new();
        for gate in self.gates.iter() {
            outputs.push(self.witness[gate.output_index]);
        }

        outputs
    }

    pub fn get_left_selectors(&self) -> Vec<u32> {
        self.gates
            .iter()
            .map(|gate| {
                if let Gate {
                    gate_type: GateType::Addition,
                    ..
                } = gate
                {
                    1
                } else {
                    0
                }
            })
            .collect()
    }

    pub fn get_right_selectors(&self) -> Vec<u32> {
        self.gates
            .iter()
            .map(|gate| {
                if let Gate {
                    gate_type: GateType::Addition,
                    ..
                } = gate
                {
                    1
                } else {
                    0
                }
            })
            .collect()
    }

    pub fn get_output_selectors(&self) -> Vec<u32> {
        self.gates
            .iter()
            .map(|_| self.field.additive_inverse(1))
            .collect()
    }

    pub fn get_multiply_selectors(&self) -> Vec<u32> {
        self.gates
            .iter()
            .map(|gate| {
                if let Gate {
                    gate_type: GateType::Multiplication,
                    ..
                } = gate
                {
                    1
                } else {
                    0
                }
            })
            .collect()
    }

    pub fn get_left_copy_constraints(&self) -> Vec<u32> {
        let mut result = Vec::new();
        for i in 0..self.gates.len() {
            let left_index = self.gates[i].left_index;

            for k in 0..self.gates.len() {
                if self.gates[k].right_index == left_index {
                    result.push(self.coset1[k]);
                } else if self.gates[k].output_index == left_index {
                    result.push(self.coset2[k]);
                } else if self.gates[k].left_index == left_index {
                    if k != i {
                        result.push(self.roots[k]);
                    }
                }
            }
        }
        result
    }

    pub fn get_right_copy_constraints(&self) -> Vec<u32> {
        let mut result = Vec::new();
        for i in 0..self.gates.len() {
            let right_index = self.gates[i].right_index;
            for k in 0..self.gates.len() {
                if self.gates[k].left_index == right_index {
                    result.push(self.roots[k]);
                } else if self.gates[k].right_index == right_index {
                    if k != i {
                        result.push(self.coset1[k]);
                    }
                } else if self.gates[k].output_index == right_index {
                    result.push(self.coset2[k]);
                }
            }
        }
        result
    }

    pub fn get_output_copy_constraints(&self) -> Vec<u32> {
        let mut result = Vec::new();
        for i in 0..self.gates.len() {
            let output_index = self.gates[i].output_index;
            for k in 0..self.gates.len() {
                if self.gates[k].left_index == output_index {
                    result.push(self.roots[k]);
                } else if self.gates[k].right_index == output_index {
                    result.push(self.coset1[k]);
                } else if self.gates[k].output_index == output_index {
                    if k != i {
                        result.push(self.coset2[k]);
                    }
                }
            }
        }
        result
    }

    pub fn get_accumulator_values(&self, degree: u32, beta: u32, gamma: u32) -> Vec<u32> {
        let mut res = Vec::new();
        let a = self.get_left_inputs();
        let b = self.get_right_inputs();
        let c = self.get_outputs();

        res.push(1);
        for round in 0..degree {
            let mut new_acc_numerator = self.field.add(
                a[round as usize],
                self.field
                    .add(self.field.multiply(beta, self.roots[round as usize]), gamma),
            );
            new_acc_numerator = self.field.multiply(
                new_acc_numerator,
                self.field.add(
                    b[round as usize],
                    self.field.add(
                        self.field.multiply(beta, self.coset1[round as usize]),
                        gamma,
                    ),
                ),
            );
            new_acc_numerator = self.field.multiply(
                new_acc_numerator,
                self.field.add(
                    c[round as usize],
                    self.field.add(
                        self.field.multiply(beta, self.coset2[round as usize]),
                        gamma,
                    ),
                ),
            );

            let mut new_acc_denominator = self.field.add(
                a[round as usize],
                self.field.add(
                    self.field.multiply(
                        beta,
                        self.circuit_polys
                            .left_copy
                            .eval(self.roots[round as usize]),
                    ),
                    gamma,
                ),
            );
            new_acc_denominator = self.field.multiply(
                new_acc_denominator,
                self.field.add(
                    b[round as usize],
                    self.field.add(
                        self.field.multiply(
                            beta,
                            self.circuit_polys
                                .right_copy
                                .eval(self.roots[round as usize]),
                        ),
                        gamma,
                    ),
                ),
            );
            new_acc_denominator = self.field.multiply(
                new_acc_denominator,
                self.field.add(
                    c[round as usize],
                    self.field.add(
                        self.field.multiply(
                            beta,
                            self.circuit_polys
                                .output_copy
                                .eval(self.roots[round as usize]),
                        ),
                        gamma,
                    ),
                ),
            );

            let acc_i = self.field.divide(new_acc_numerator, new_acc_denominator);
            res.push(self.field.multiply(res[(round) as usize], acc_i));
        }

        res
    }
}

pub enum GateType {
    Addition,
    Multiplication,
}

pub struct Gate {
    pub left_index: usize,
    pub right_index: usize,
    pub output_index: usize,
    pub gate_type: GateType,
}

#[derive(Default)]
pub struct CircuitPolynomials {
    pub left: Polynomial,
    pub right: Polynomial,
    pub output: Polynomial,
    pub left_selector: Polynomial,
    pub right_selector: Polynomial,
    pub output_selector: Polynomial,
    pub multiply_selector: Polynomial,
    pub left_copy: Polynomial,
    pub right_copy: Polynomial,
    pub output_copy: Polynomial,
    pub z_h: Polynomial,
    pub acc: Polynomial,
}
