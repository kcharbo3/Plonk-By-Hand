use crate::field::Field;
use crate::math::ecc::{CurvePoint, ExtensionCurvePoint, ECC};

#[derive(Clone)]
pub struct SRS {
    pub g_1: CurvePoint,
    pub g_2: ExtensionCurvePoint,
    pub g_1_points: Vec<CurvePoint>,
    pub g_2_points: Vec<ExtensionCurvePoint>,
    pub degree: u32,
    pub s: u32,
    pub scalar_field: Field,
    pub ecc: ECC,
}

impl SRS {
    pub fn generate_g_1_points(&mut self) {
        let mut g_1_points = Vec::new();
        for i in 0..(self.degree + 3) {
            g_1_points.push(
                self.ecc
                    .multiply(self.scalar_field.exponent(self.s, i), &self.g_1),
            );
        }

        self.g_1_points = g_1_points;
    }

    pub fn generate_g_2_points(&mut self) {
        self.g_2_points.push(self.g_2.clone());
        self.g_2_points.push(self.ecc.double_extension(&self.g_2));
    }

    pub fn copy(&self) -> SRS {
        SRS {
            g_1: self.g_1.clone(),
            g_2: self.g_2.clone(),
            g_1_points: self.g_1_points.clone(),
            g_2_points: self.g_2_points.clone(),
            degree: self.degree,
            s: self.s,
            scalar_field: self.scalar_field.clone(),
            ecc: ECC {
                field: self.ecc.field.clone(),
            },
        }
    }
}

#[test]
fn test_g_1_points() {
    let mut srs = SRS {
        g_1: CurvePoint::new(1, 2),
        g_2: ExtensionCurvePoint::new(31, 36, true),
        g_1_points: Vec::new(),
        g_2_points: Vec::new(),
        degree: 4,
        s: 2,
        scalar_field: Field { order: 17 },
        ecc: ECC {
            field: Field { order: 101 },
        },
    };

    srs.generate_g_1_points();
    let expected_points = vec![
        CurvePoint { x: 1, y: 2 },
        CurvePoint { x: 68, y: 74 },
        CurvePoint { x: 65, y: 98 },
        CurvePoint { x: 18, y: 49 },
        CurvePoint { x: 1, y: 99 },
        CurvePoint { x: 68, y: 27 },
        CurvePoint { x: 65, y: 3 },
    ];

    for i in 0..srs.g_1_points.len() {
        assert!(CurvePoint::equals(&srs.g_1_points[i], &expected_points[i]));
    }
}

#[test]
fn test_g_2_points() {
    let mut srs = SRS {
        g_1: CurvePoint::new(1, 2),
        g_2: ExtensionCurvePoint::new(36, 31, true),
        g_1_points: Vec::new(),
        g_2_points: Vec::new(),
        degree: 4,
        s: 2,
        scalar_field: Field { order: 17 },
        ecc: ECC {
            field: Field { order: 101 },
        },
    };

    srs.generate_g_2_points();

    let expected_first_point = ExtensionCurvePoint::new(36, 31, true);
    let expected_second_point = ExtensionCurvePoint::new(90, 82, true);
    assert!(ExtensionCurvePoint::equals(
        &srs.g_2_points[0],
        &expected_first_point
    ));
    assert!(ExtensionCurvePoint::equals(
        &srs.g_2_points[1],
        &expected_second_point
    ));
}
