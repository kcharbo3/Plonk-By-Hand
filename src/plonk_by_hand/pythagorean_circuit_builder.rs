use crate::math::field;
use crate::math::roots_of_unity;
use crate::plonk_by_hand::circuit;

pub struct PythagoreanCircuit {
    pub circuit: circuit::Circuit,
}

impl PythagoreanCircuit {
    pub fn new(field: field::Field) -> PythagoreanCircuit {
        let roots = roots_of_unity::get_roots_of_unity(4, &field);

        PythagoreanCircuit {
            circuit: circuit::Circuit {
                gates: Vec::new(),
                inputs: Vec::new(),
                witness: vec![0; 6],
                field: field.clone(),
                coset1: roots_of_unity::get_coset(2, &roots, &field),
                coset2: roots_of_unity::get_coset(3, &roots, &field),
                roots,
                circuit_polys: Default::default(),
            },
        }
    }

    pub fn build_circuit_with_inputs(&mut self, inputs: Vec<u32>) {
        // inputs[0]: X_0, inputs[1]: X_2, inputs[2]: X_4
        self.circuit
            .set_inputs(vec![(inputs[0], 0), (inputs[1], 2), (inputs[2], 4)]);

        self.build_circuit();

        self.circuit.compute_witness();

        // private polynomials
        self.circuit.build_polynomials_with_input();
    }

    pub fn build_circuit(&mut self) {
        // X_0 * X_0 = X_1
        self.circuit.insert_gate(circuit::Gate {
            left_index: 0,
            right_index: 0,
            output_index: 1,
            gate_type: circuit::GateType::Multiplication,
        });

        // X_2 * X_2 = X_3
        self.circuit.insert_gate(circuit::Gate {
            left_index: 2,
            right_index: 2,
            output_index: 3,
            gate_type: circuit::GateType::Multiplication,
        });

        // X_4 * X_4 = X_5
        self.circuit.insert_gate(circuit::Gate {
            left_index: 4,
            right_index: 4,
            output_index: 5,
            gate_type: circuit::GateType::Multiplication,
        });

        // X_1 + X_3 = X_5
        self.circuit.insert_gate(circuit::Gate {
            left_index: 1,
            right_index: 3,
            output_index: 5,
            gate_type: circuit::GateType::Addition,
        });

        // public polynomials
        self.circuit.build_polynomials();
    }

    pub fn build_acc(&mut self, beta: u32, gamma: u32) {
        self.circuit.build_accumulator(beta, gamma);
    }
}

#[test]
fn test_circuit() {
    let field = field::Field { order: 17 };

    let mut pythagorean_circuit = PythagoreanCircuit::new(field);
    pythagorean_circuit.build_circuit_with_inputs(vec![3, 4, 5]);
    assert_eq!(
        vec![3, 4, 5, 9],
        pythagorean_circuit.circuit.get_left_inputs()
    );
    assert_eq!(
        vec![3, 4, 5, 16],
        pythagorean_circuit.circuit.get_right_inputs()
    );
    assert_eq!(vec![9, 16, 8, 8], pythagorean_circuit.circuit.get_outputs());

    assert_eq!(
        vec![0, 0, 0, 1],
        pythagorean_circuit.circuit.get_left_selectors()
    );
    assert_eq!(
        vec![0, 0, 0, 1],
        pythagorean_circuit.circuit.get_right_selectors()
    );
    assert_eq!(
        vec![16, 16, 16, 16],
        pythagorean_circuit.circuit.get_output_selectors()
    );
    assert_eq!(
        vec![1, 1, 1, 0],
        pythagorean_circuit.circuit.get_multiply_selectors()
    );

    assert_eq!(
        vec![2, 8, 15, 3],
        pythagorean_circuit.circuit.get_left_copy_constraints()
    );
    assert_eq!(
        vec![1, 4, 16, 12],
        pythagorean_circuit.circuit.get_right_copy_constraints()
    );
    assert_eq!(
        vec![13, 9, 5, 14],
        pythagorean_circuit.circuit.get_output_copy_constraints()
    );

    let left_poly = pythagorean_circuit.circuit.circuit_polys.left;
    assert_eq!(vec![1, 13, 3, 3], left_poly.coefficients);
    let right_poly = pythagorean_circuit.circuit.circuit_polys.right;
    assert_eq!(vec![7, 3, 14, 13], right_poly.coefficients);
    let output_poly = pythagorean_circuit.circuit.circuit_polys.output;
    assert_eq!(vec![6, 5, 11, 4], output_poly.coefficients);

    let left_selector_poly = pythagorean_circuit.circuit.circuit_polys.left_selector;
    assert_eq!(vec![13, 1, 4, 16], left_selector_poly.coefficients);
    let right_selector_poly = pythagorean_circuit.circuit.circuit_polys.right_selector;
    assert_eq!(vec![13, 1, 4, 16], right_selector_poly.coefficients);
    let output_selector_poly = pythagorean_circuit.circuit.circuit_polys.output_selector;
    assert_eq!(vec![16, 0, 0, 0], output_selector_poly.coefficients);
    let multiply_selector_poly = pythagorean_circuit.circuit.circuit_polys.multiply_selector;
    assert_eq!(vec![5, 16, 13, 1], multiply_selector_poly.coefficients);

    let left_copy_poly = pythagorean_circuit.circuit.circuit_polys.left_copy;
    assert_eq!(vec![7, 13, 10, 6], left_copy_poly.coefficients);
    let right_copy_poly = pythagorean_circuit.circuit.circuit_polys.right_copy;
    assert_eq!(vec![4, 0, 13, 1], right_copy_poly.coefficients);
    let output_copy_poly = pythagorean_circuit.circuit.circuit_polys.output_copy;
    assert_eq!(vec![6, 7, 3, 14], output_copy_poly.coefficients);

    assert_eq!(
        vec![16, 0, 0, 0, 1],
        pythagorean_circuit.circuit.circuit_polys.z_h.coefficients
    );
}
