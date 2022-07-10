use crate::Field;

#[derive(Debug)]
pub struct ComplexScalar {
    pub constant: u32,
    pub u_term: u32,
}

impl ComplexScalar {
    pub fn multiply(field: &Field, c1: &ComplexScalar, c2: &ComplexScalar) -> ComplexScalar {
        let mut constant = field.multiply(c1.constant, c2.constant);
        constant = field.add(
            constant,
            field.multiply(
                field.additive_inverse(2),
                field.multiply(c1.u_term, c2.u_term),
            ),
        );
        let u_term = field.add(
            field.multiply(c1.constant, c2.u_term),
            field.multiply(c1.u_term, c2.constant),
        );

        ComplexScalar { constant, u_term }
    }

    pub fn add(field: &Field, c1: &ComplexScalar, c2: &ComplexScalar) -> ComplexScalar {
        ComplexScalar {
            constant: field.add(c1.constant, c2.constant),
            u_term: field.add(c1.u_term, c2.u_term),
        }
    }

    pub fn equals(c1: &ComplexScalar, c2: &ComplexScalar) -> bool {
        c1.constant == c2.constant && c1.u_term == c2.u_term
    }

    pub fn exponent(field: &Field, c: &ComplexScalar, exponent: u32) -> ComplexScalar {
        let mut curr = 2;
        let mut res = ComplexScalar {
            constant: c.constant,
            u_term: c.u_term,
        };
        // 1 -> 2 -> 4 -> 8 -> 16 -> 32 -> 64 -> 128 -> 256 -> 512
        while curr < exponent {
            res = ComplexScalar::multiply(field, &res, &res);
            curr *= 2;
        }

        curr /= 2;
        while curr < exponent {
            res = ComplexScalar::multiply(field, c, &res);
            curr += 1;
        }

        res
    }
}
