use crate::math::field::Field;
use crate::math::matrix;

pub struct Polynomial {
    pub degree: u32,
    pub coefficients: Vec<u32>,
    pub field: Field,
}

#[derive(Debug)]
pub struct Point {
    pub x: u32,
    pub y: u32,
}

impl Default for Polynomial {
    fn default() -> Polynomial {
        Polynomial {
            degree: 0,
            coefficients: Vec::new(),
            field: Field { order: 0 },
        }
    }
}

impl Polynomial {
    pub fn eval(&self, x: u32) -> u32 {
        let mut res = self.coefficients[0];
        for i in 1..(self.degree + 1) {
            let x_i = self.field.exponent(x, i);
            res = self
                .field
                .add(res, self.field.multiply(self.coefficients[i as usize], x_i));
        }

        res
    }

    pub fn scalar_multiply(&self, scalar: u32) -> Polynomial {
        let new_coefficients = self
            .coefficients
            .iter()
            .map(|coefficient| self.field.multiply(scalar, *coefficient))
            .collect();

        Polynomial::remove_zeros(&Polynomial {
            degree: self.degree,
            coefficients: new_coefficients,
            field: self.field.clone(),
        })
    }

    // 0*0, 1*0 0*1, 2*0 0*2 1*1, ...
    // TODO: impl FFT
    pub fn poly_multiply(poly1: &Polynomial, poly2: &Polynomial) -> Polynomial {
        let new_degree = poly1.degree + poly2.degree;
        let mut new_coefficients = vec![0; (new_degree + 1) as usize];

        for i in 0..=new_degree {
            for k in 0..=poly1.degree {
                for j in 0..=poly2.degree {
                    if k + j == i {
                        new_coefficients[i as usize] = poly1.field.add(
                            new_coefficients[i as usize],
                            poly1.field.multiply(
                                poly1.coefficients[k as usize],
                                poly2.coefficients[j as usize],
                            ),
                        );
                    }
                }
            }
        }

        Polynomial::remove_zeros(&Polynomial {
            degree: new_degree,
            coefficients: new_coefficients,
            field: poly1.field.clone(),
        })
    }

    pub fn poly_add(poly1: &Polynomial, poly2: &Polynomial) -> Polynomial {
        let mut large = poly1;
        let mut small = poly2;
        if poly1.degree < poly2.degree {
            large = poly2;
            small = poly1;
        }

        let mut res = large.copy();
        for i in 0..small.coefficients.len() {
            res.coefficients[i] = poly1
                .field
                .add(large.coefficients[i], small.coefficients[i]);
        }

        Polynomial::remove_zeros(&res)
    }

    pub fn poly_subtract(poly1: &Polynomial, poly2: &Polynomial) -> Polynomial {
        Polynomial::poly_add(poly1, &Polynomial::additive_inverse_poly(poly2))
    }

    pub fn additive_inverse_poly(poly: &Polynomial) -> Polynomial {
        let mut new_coefficients = vec![0; (poly.degree + 1) as usize];
        for i in 0..=poly.degree {
            new_coefficients[i as usize] =
                poly.field.additive_inverse(poly.coefficients[i as usize]);
        }

        Polynomial {
            degree: poly.degree,
            coefficients: new_coefficients,
            field: poly.field.clone(),
        }
    }

    pub fn poly_divide(poly1: &Polynomial, poly2: &Polynomial) -> (Polynomial, Polynomial) {
        //function n / d:
        // 2      require d ≠ 0
        // 3      (q, r) ← (0, n)            # At each step n = d × q + r
        // 4      while r ≠ 0 AND degree(r) ≥ degree(d):
        // 5         t ← lead(r)/lead(d)     # Divide the leading terms
        // 6         (q, r) ← (q + t, r - (t * d))
        // 7      return (q, r)

        let mut q = Polynomial {
            field: poly1.field.clone(),
            coefficients: Vec::new(),
            degree: 0,
        };

        let mut r = Polynomial {
            field: poly1.field.clone(),
            coefficients: poly1.coefficients.clone(),
            degree: poly1.degree,
        };

        while (r.coefficients[0] != 0 && r.coefficients.len() != 1) && r.degree >= poly2.degree {
            let new_coefficient = poly1.field.divide(
                r.coefficients[r.coefficients.len() - 1],
                poly2.coefficients[poly2.coefficients.len() - 1],
            );

            let new_degree = r.degree - poly2.degree;
            let mut new_coefficient_vec = vec![0; (new_degree + 1) as usize];
            new_coefficient_vec[new_degree as usize] = new_coefficient;

            let t = Polynomial {
                degree: new_degree,
                coefficients: new_coefficient_vec,
                field: poly1.field.clone(),
            };

            q = Polynomial::poly_add(&q, &t);
            r = Polynomial::poly_add(
                &r,
                &Polynomial::additive_inverse_poly(&Polynomial::poly_multiply(&t, poly2)),
            );
        }

        (q, r)
    }

    pub fn from_scalar(scalar: u32, field: Field) -> Polynomial {
        Polynomial {
            degree: 0,
            coefficients: vec![scalar],
            field,
        }
    }

    fn remove_zeros(poly: &Polynomial) -> Polynomial {
        let mut first_non_zero_index = 0;
        for i in (0..poly.coefficients.len()).rev() {
            if poly.coefficients[i] != 0 {
                first_non_zero_index = i;
                break;
            }
        }

        Polynomial {
            degree: first_non_zero_index as u32,
            coefficients: poly.coefficients[0..=first_non_zero_index].to_vec(),
            field: poly.field.clone(),
        }
    }

    pub fn copy(&self) -> Polynomial {
        Polynomial {
            degree: self.degree,
            coefficients: self.coefficients.clone(),
            field: self.field.clone(),
        }
    }

    pub fn create_poly_from_points(points: Vec<Point>, field: &Field) -> Polynomial {
        // todo: compute inverse matrix...
        let mut y_vals = Vec::new();
        let mut x_vals = Vec::new();
        for point in points.iter() {
            x_vals.push(point.x);
            y_vals.push(point.y);
        }

        let interpolation_matrix = Polynomial::get_interpolation_matrix(x_vals, field);
        let inv_matrix = matrix::get_inverse_matrix_4x4(&interpolation_matrix, field);
        let coefficients = matrix::matrix_multiply_4x4_1x4(&inv_matrix, &y_vals, field);

        Polynomial {
            degree: (points.len() - 1) as u32,
            coefficients,
            field: field.clone(),
        }
    }

    fn get_interpolation_matrix(x_vals: Vec<u32>, field: &Field) -> Vec<Vec<u32>> {
        let mut res = Vec::new();
        for i in 0..x_vals.len() {
            res.push(Polynomial::get_interpolation_row(
                x_vals[i],
                x_vals.len() as u32,
                field,
            ));
        }

        res
    }

    fn get_interpolation_row(x: u32, degree: u32, field: &Field) -> Vec<u32> {
        let mut res = Vec::new();
        let mut curr;
        for i in 0..(degree + 1) {
            curr = field.exponent(x, i);
            res.push(curr);
        }

        res
    }

    pub fn get_lagrange_1_poly(field: &Field, degree: u32, roots: &Vec<u32>) -> Polynomial {
        let mut points = Vec::new();
        for i in 0..=degree {
            let y = if i == 0 { 1 } else { 0 };
            points.push(Point {
                x: roots[i as usize],
                y,
            });
        }

        Polynomial::create_poly_from_points(points, field)
    }
}

#[test]
fn test_create_poly_from_points() {
    let field = Field { order: 17 };

    let points = vec![
        Point { x: 1, y: 3 },
        Point { x: 4, y: 4 },
        Point { x: 16, y: 5 },
        Point { x: 13, y: 9 },
    ];

    let new_poly = Polynomial::create_poly_from_points(points, &field);
    let expected_coefficients: Vec<u32> = vec![1, 13, 3, 3];

    for i in 0..new_poly.coefficients.len() {
        assert_eq!(new_poly.coefficients[i], expected_coefficients[i]);
    }

    assert_eq!(new_poly.degree, 3);
}

#[test]
fn test_poly_multiply() {
    let field = Field { order: 17 };
    let poly1 = Polynomial {
        degree: 2,
        coefficients: vec![16, 2, 1],
        field: field.clone(),
    };

    let poly2 = Polynomial {
        degree: 2,
        coefficients: vec![6, 14, 2],
        field,
    };

    let actual_poly = Polynomial::poly_multiply(&poly1, &poly2);
    let expected_coefficients = vec![11, 15, 15, 1, 2];
    assert_eq!(actual_poly.coefficients, expected_coefficients);
}

#[test]
fn test_poly_divide() {
    let field = Field { order: 17 };

    let poly1 = Polynomial {
        degree: 2,
        coefficients: vec![2, 1, 2],
        field: field.clone(),
    };

    let poly2 = Polynomial {
        degree: 1,
        coefficients: vec![3, 1],
        field,
    };

    let (quotient, remainder) = Polynomial::poly_divide(&poly1, &poly2);
    let expected_quotient = vec![12, 2];
    let expected_remainder = vec![0];
    assert_eq!(quotient.coefficients, expected_quotient);
    assert_eq!(remainder.coefficients, expected_remainder);
}
