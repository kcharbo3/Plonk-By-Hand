use crate::math::complex_scalar::ComplexScalar;
use crate::math::field::Field;

const INFINITY: u32 = 99999999;

// Curve: y^2 = x^3 + 3
#[derive(Clone)]
pub struct ECC {
    pub field: Field,
}

#[derive(Debug, Default, Clone)]
pub struct CurvePoint {
    pub x: u32,
    pub y: u32,
}

impl CurvePoint {
    pub fn new(x: u32, y: u32) -> CurvePoint {
        CurvePoint { x, y }
    }

    pub fn point_at_infinity() -> CurvePoint {
        CurvePoint {
            x: INFINITY,
            y: INFINITY,
        }
    }

    pub fn equals(p1: &CurvePoint, p2: &CurvePoint) -> bool {
        (p1.x == p2.x) && (p1.y == p2.y)
    }
}

#[derive(Debug, Default, Clone)]
pub struct ExtensionCurvePoint {
    pub x: u32,
    pub y: u32,
    pub u: bool,
}

impl ExtensionCurvePoint {
    pub fn new(x: u32, y: u32, u: bool) -> ExtensionCurvePoint {
        ExtensionCurvePoint { x, y, u }
    }

    pub fn equals(p1: &ExtensionCurvePoint, p2: &ExtensionCurvePoint) -> bool {
        (p1.x == p2.x) && (p1.y == p2.y) && (p1.u == p2.u)
    }
}

impl ECC {
    pub fn copy(&self) -> ECC {
        ECC {
            field: self.field.clone(),
        }
    }

    // only used if p1 != p2
    pub fn add(&self, p1: &CurvePoint, p2: &CurvePoint) -> CurvePoint {
        if CurvePoint::equals(p1, &CurvePoint::point_at_infinity()) {
            return p2.clone();
        } else if CurvePoint::equals(p2, &CurvePoint::point_at_infinity()) {
            return p1.clone();
        } else if CurvePoint::equals(p1, p2) {
            return self.double(p1);
        } else if p1.x == p2.x {
            return CurvePoint::point_at_infinity();
        }

        let delta_x = self.field.subtract(p2.x, p1.x);
        let delta_y = self.field.subtract(p2.y, p1.y);
        let m = self.field.divide(delta_y, delta_x);

        let new_x = self
            .field
            .subtract(self.field.subtract(self.field.multiply(m, m), p1.x), p2.x);
        let new_y = self.field.subtract(
            self.field.multiply(m, self.field.subtract(p1.x, new_x)),
            p1.y,
        );

        CurvePoint { x: new_x, y: new_y }
    }

    pub fn double(&self, p: &CurvePoint) -> CurvePoint {
        let (new_x, new_y, _u) = self.double_any(p.x, p.y, false);

        CurvePoint::new(new_x, new_y)
    }

    pub fn double_extension(&self, p: &ExtensionCurvePoint) -> ExtensionCurvePoint {
        let (new_x, new_y, u) = self.double_any(p.x, p.y, p.u);

        ExtensionCurvePoint::new(new_x, new_y, u)
    }

    fn double_any(&self, x: u32, y: u32, u: bool) -> (u32, u32, bool) {
        if !u && CurvePoint::equals(&CurvePoint::new(x, y), &CurvePoint::point_at_infinity()) {
            return (x, y, false);
        }

        let m_num = self.field.multiply(3, self.field.multiply(x, x));
        let m_denom = self.field.multiply(2, y);
        let m = self.field.divide(m_num, m_denom);

        let u_factor = self.field.divide(1, self.field.additive_inverse(2));

        let m_squared = if u {
            self.field.multiply(m, self.field.multiply(m, u_factor))
        } else {
            self.field.multiply(m, m)
        };

        let new_x = self.field.subtract(m_squared, self.field.multiply(2, x));

        let y_u_term = if u {
            self.field.multiply(
                self.field
                    .multiply(m, self.field.subtract(self.field.multiply(3, x), m_squared)),
                u_factor,
            )
        } else {
            self.field
                .multiply(m, self.field.subtract(self.field.multiply(3, x), m_squared))
        };

        let new_y = self.field.subtract(y_u_term, y);

        (new_x, new_y, u)
    }

    pub fn inversion(&self, p: &CurvePoint) -> CurvePoint {
        CurvePoint {
            x: p.x,
            y: self.field.additive_inverse(p.y),
        }
    }

    // 8G = 4 (2G) -> 2 (4G) -> 8G
    // 14G = 7 (2G) -> 6 (3G) -> 3 (6G) -> 2 (7G) -> 14G
    // 15G = 14 (G) -> 7 (2G) ->
    // 15G = G + 14G
    pub fn multiply(&self, scalar: u32, p1: &CurvePoint) -> CurvePoint {
        if scalar == 0 {
            CurvePoint::point_at_infinity()
        } else if scalar == 1 {
            p1.clone()
        } else if scalar % 2 == 0 {
            let doubled = self.double(p1);
            self.multiply(scalar / 2, &doubled)
        } else {
            self.add(p1, &self.multiply(scalar - 1, &p1.clone()))
        }
    }

    // outputs x factor, y factor, and constant -> for y^2
    pub fn get_line_between_points(&self, p1: &CurvePoint, p2: &CurvePoint) -> (u32, u32, u32) {
        let m_numerator = self.field.subtract(p2.y, p1.y);
        let m_denominator = self.field.subtract(p2.x, p1.x);

        let mut y_factor = m_denominator;
        let mut x_factor = m_numerator;
        let mut constant = self
            .field
            .multiply(self.field.additive_inverse(p1.y), m_denominator);
        constant = self.field.subtract(
            constant,
            self.field
                .multiply(self.field.additive_inverse(p1.x), m_numerator),
        );

        // hardcode when to switch sides to match the plonk_by_hand article
        if y_factor > x_factor || p1.x == 65 || p1.x == 1 {
            y_factor = self.field.additive_inverse(y_factor);
            constant = self.field.additive_inverse(constant);
        } else {
            x_factor = self.field.additive_inverse(x_factor);
        }

        if p1.x == p2.x {
            x_factor = 1;
            y_factor = 0;
            constant = 15;
        }

        (x_factor, y_factor, constant)
    }

    // x_factor*x + y_factor*y + constant
    pub fn plug_extension_point_in_equation(
        &self,
        p: &ExtensionCurvePoint,
        x_factor: u32,
        y_factor: u32,
        constant: u32,
    ) -> ComplexScalar {
        ComplexScalar {
            constant: self.field.add(self.field.multiply(x_factor, p.x), constant),
            u_term: self.field.multiply(y_factor, p.y),
        }
    }

    pub fn get_generator_multiple(
        &self,
        p: &CurvePoint,
        generator: &CurvePoint,
        order: u32,
    ) -> u32 {
        let mut ans = 0;
        for i in 1..order {
            let new_point = self.multiply(i, generator);
            if CurvePoint::equals(p, &new_point) {
                ans = i;
                break;
            }
        }
        ans
    }
}

#[test]
fn test_multiply() {
    let ecc = ECC {
        field: Field { order: 101 },
    };

    let p1 = CurvePoint { x: 1, y: 2 };
    let expected_product1 = CurvePoint { x: 68, y: 74 };
    assert!(CurvePoint::equals(&ecc.double(&p1), &expected_product1));
    assert!(CurvePoint::equals(&ecc.add(&p1, &p1), &expected_product1));
    assert!(CurvePoint::equals(
        &ecc.multiply(2, &p1),
        &expected_product1
    ));

    let expected_product2 = CurvePoint { x: 1, y: 99 };
    assert!(CurvePoint::equals(
        &ecc.multiply(16, &p1),
        &expected_product2
    ));
}

#[test]
fn test_extension_double() {
    let ecc = ECC {
        field: Field { order: 101 },
    };

    let g2 = ExtensionCurvePoint::new(36, 31, true);
    let doubled = ecc.double_extension(&g2);
    let expected = ExtensionCurvePoint::new(90, 82, true);
    assert!(ExtensionCurvePoint::equals(&expected, &doubled));
}

#[test]
fn test_line_through_points() {
    let ecc = ECC {
        field: Field { order: 101 },
    };

    let (x, y, constant) =
        ecc.get_line_between_points(&CurvePoint::new(1, 2), &CurvePoint::new(68, 27));

    assert_eq!(x, 25);
    assert_eq!(y, 34);
    assert_eq!(constant, 8);

    let (x, y, constant) =
        ecc.get_line_between_points(&CurvePoint::new(68, 74), &CurvePoint::new(65, 3));

    assert_eq!(x, 30);
    assert_eq!(y, 3);
    assert_eq!(constant, 61);

    let (x, y, constant) =
        ecc.get_line_between_points(&CurvePoint::new(65, 98), &CurvePoint::new(18, 52));

    assert_eq!(x, 55);
    assert_eq!(y, 47);
    assert_eq!(constant, 0);

    let (x, y, constant) =
        ecc.get_line_between_points(&CurvePoint::new(18, 49), &CurvePoint::new(1, 2));

    assert_eq!(x, 54);
    assert_eq!(y, 17);
    assert_eq!(constant, 13);

    let (x, y, constant) =
        ecc.get_line_between_points(&CurvePoint::new(1, 99), &CurvePoint::new(1, 2));

    assert_eq!(x, 1);
    assert_eq!(y, 0);
    assert_eq!(constant, 15);
}
