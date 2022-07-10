use crate::reed_solomon::ReedSolomon;
use crate::field::Field;

use ndarray::{array, Array1, Array2, ArrayView1, ArrayView2};

pub struct Freivald {
    pub reed: ReedSolomon
}

impl Freivald {
    pub fn verify(&self, a: &ArrayView2::<u32>, b: &ArrayView2::<u32>, c: &ArrayView2::<u32>, challenge: u32) -> bool {

        let y = self.compute_y(c, challenge);
        let z = self.compute_z(a, b, challenge);

        let mut result = true;
        for (index, val) in y.indexed_iter() {
            if y[index] != z[index] {
                result = false;
                break;
            }
        }

        result
    }

    // y = Cx
    pub fn compute_y(&self, c: &ArrayView2::<u32>, challenge: u32) -> Array1::<u32> {

        let mut result: Vec<u32> = Vec::new();

        for row in c.rows() {
            result.push(self.reed.finger_print(&row.view() as &ArrayView1<u32>, challenge));
        }

        Array1::<u32>::from_vec(result)
    }

    // z = A * Bx
    pub fn compute_z(&self, a: &ArrayView2::<u32>, b: &ArrayView2::<u32>, challenge: u32) -> Array1::<u32> {
        let mut b_mul_x: Vec<u32> = Vec::new();
        for row in b.rows() {
            b_mul_x.push(self.reed.finger_print(&row.view() as &ArrayView1<u32>, challenge));
        }

        let b_mul_x_arr: Array1<u32> = Array1::<u32>::from_vec(b_mul_x);

        let mut result : Vec<u32> = Vec::new();
        for row in a.rows() {
            &result.push(self.reed.field.vector_dot(
                &row.view() as &ArrayView1<u32>,
                &b_mul_x_arr.view()));
        }

        Array1::<u32>::from_vec(result)
    }

}

#[test]
fn test_verify() {
    let freivald = Freivald {
        reed: ReedSolomon {
            field: Field {
                order: 9999999
            }
        }
    };

    let a: Array2::<u32> = array![[1, 2], [3, 4]];
    let b: Array2::<u32> = array![[2, 3], [4, 5]];
    let c = &a.dot(&b);

    assert_eq!(freivald.verify(&a.view(), &b.view(), &c.view(), 13), true);
}