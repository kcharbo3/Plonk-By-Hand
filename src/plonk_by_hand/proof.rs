use crate::math::ecc::CurvePoint;
use crate::math::polynomial::Polynomial;

#[derive(Default)]
pub struct Proof {
    pub a: CurvePoint,
    pub b: CurvePoint,
    pub c: CurvePoint,
    pub z: CurvePoint,
    pub t_lo: CurvePoint,
    pub t_mid: CurvePoint,
    pub t_hi: CurvePoint,
    pub w: CurvePoint,
    pub wz: CurvePoint,
    pub a_bar: u32,
    pub b_bar: u32,
    pub c_bar: u32,
    pub left_copy_bar: u32,
    pub right_copy_bar: u32,
    pub r_bar: u32,
    pub z_bar: u32,
}

#[derive(Default)]
pub struct ProverPolys {
    pub a: Polynomial,
    pub b: Polynomial,
    pub c: Polynomial,
    pub z: Polynomial,
    pub t: Polynomial,
    pub t_lo: Polynomial,
    pub t_mid: Polynomial,
    pub t_hi: Polynomial,
    pub r: Polynomial,
    pub w: Polynomial,
    pub wz: Polynomial,
}

#[derive(Default, Clone)]
pub struct OpeningEvals {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub left_copy: u32,
    pub right_copy: u32,
    pub t: u32,
    pub z: u32,
    pub r: u32,
}
