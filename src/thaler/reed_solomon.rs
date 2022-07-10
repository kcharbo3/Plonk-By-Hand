use crate::field::Field;
use ndarray::{Array1, ArrayView1};

pub struct ReedSolomon {
    pub field: Field
}

impl ReedSolomon {
    pub fn finger_print(&self, values: &ArrayView1<u32>, challenge: u32) -> u32 {
        let challenge_array = self.get_challenge_array(challenge, values.len() as u32);

        self.field.vector_dot(values, &challenge_array.view())
    }

    pub fn get_challenge_array(&self, challenge: u32, degree: u32) -> Array1::<u32> {
        let mut challenge_array: Vec<u32> = Vec::new();
        let mut val = 1;
        for index in 0..degree {
            challenge_array.push(val);
            val = self.field.multiply(val, challenge);
        }

        Array1::<u32>::from_vec(challenge_array)
    }
}