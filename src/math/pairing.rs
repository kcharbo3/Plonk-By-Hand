use crate::math::ecc::ComplexScalar;
use crate::plonk_by_hand::constants;
use crate::{CurvePoint, ExtensionCurvePoint, Field, ECC};

pub struct Pairing {
    pub r: u32,
    pub ecc: ECC,
}

impl Pairing {
    pub fn get_base_pairing(&self, q: &ExtensionCurvePoint, p: &CurvePoint) -> ComplexScalar {
        let f17 = self.get_f_17(q, p);
        let exponent = (u32::pow(self.ecc.field.order, 2) - 1) / (self.r);

        ComplexScalar::exponent(&self.ecc.field, &f17, exponent)
    }

    fn get_f_17(&self, q: &ExtensionCurvePoint, p: &CurvePoint) -> ComplexScalar {
        let p2minus = self.ecc.inversion(&self.ecc.multiply(2, p));
        let (f2_x, f2_y, f2_constant) = self.ecc.get_line_between_points(p, &p2minus);
        let f2 = self
            .ecc
            .plug_extension_point_in_equation(q, f2_x, f2_y, f2_constant);

        let p2 = self.ecc.multiply(2, p);
        let p4minus = self.ecc.inversion(&self.ecc.multiply(4, p));
        let (f4_x, f4_y, f4_constant) = self.ecc.get_line_between_points(&p2, &p4minus);
        let f4 = ComplexScalar::multiply(
            &self.ecc.field,
            &ComplexScalar::multiply(&self.ecc.field, &f2, &f2),
            &self
                .ecc
                .plug_extension_point_in_equation(q, f4_x, f4_y, f4_constant),
        );

        let p4 = self.ecc.multiply(4, p);
        let p8minus = self.ecc.inversion(&self.ecc.multiply(8, p));
        let (f8_x, f8_y, f8_constant) = self.ecc.get_line_between_points(&p4, &p8minus);
        let f8 = ComplexScalar::multiply(
            &self.ecc.field,
            &ComplexScalar::multiply(&self.ecc.field, &f4, &f4),
            &self
                .ecc
                .plug_extension_point_in_equation(q, f8_x, f8_y, f8_constant),
        );

        let p8 = self.ecc.multiply(8, p);
        let p16minus = self.ecc.inversion(&self.ecc.multiply(16, p));
        let (f16_x, f16_y, f16_constant) = self.ecc.get_line_between_points(&p8, &p16minus);
        let f16 = ComplexScalar::multiply(
            &self.ecc.field,
            &ComplexScalar::multiply(&self.ecc.field, &f8, &f8),
            &self
                .ecc
                .plug_extension_point_in_equation(q, f16_x, f16_y, f16_constant),
        );

        let p16 = self.ecc.multiply(16, p);
        let pminus = self.ecc.inversion(&p);
        let (f17_x, f17_y, f17_constant) = self.ecc.get_line_between_points(&p16, &pminus);
        let f17 = ComplexScalar::multiply(
            &self.ecc.field,
            &f16,
            &self
                .ecc
                .plug_extension_point_in_equation(q, f17_x, f17_y, f17_constant),
        );
        f17
    }
}

#[test]
fn test_f_17() {
    let pairing = Pairing {
        r: 17,
        ecc: ECC {
            field: Field { order: 101 },
        },
    };

    let mut srs = constants::SRS_BY_HAND.clone();
    srs.generate_g_1_points();
    srs.generate_g_2_points();

    let f17 = pairing.get_f_17(&srs.g_2_points[1], &srs.g_1_points[0]);
    assert!(ComplexScalar::equals(
        &f17,
        &ComplexScalar {
            constant: 68,
            u_term: 47
        }
    ));
}

#[test]
fn test_get_base_pairing() {
    let pairing = Pairing {
        r: 17,
        ecc: ECC {
            field: Field { order: 101 },
        },
    };

    let mut srs = constants::SRS_BY_HAND.clone();
    srs.generate_g_1_points();
    srs.generate_g_2_points();

    let base_pairing = pairing.get_base_pairing(&srs.g_2_points[1], &srs.g_1_points[0]);
    assert!(ComplexScalar::equals(
        &base_pairing,
        &ComplexScalar {
            constant: 97,
            u_term: 89
        }
    ));

    let second_base_pairing = pairing.get_base_pairing(&srs.g_2_points[0], &srs.g_1_points[0]);
    assert!(ComplexScalar::equals(
        &second_base_pairing,
        &ComplexScalar {
            constant: 7,
            u_term: 28
        }
    ));
}
