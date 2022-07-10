extern crate core;

pub mod math;
pub mod plonk_by_hand;

use crate::field::Field;
use crate::math::ecc::{CurvePoint, ExtensionCurvePoint, ECC};
use crate::plonk_by_hand::constants;
use crate::plonk_by_hand::prover::Prover;
use crate::plonk_by_hand::public_coin::PublicCoin;
use crate::plonk_by_hand::pythagorean_circuit_builder::PythagoreanCircuit;
use crate::plonk_by_hand::structured_reference_string::SRS;
use crate::plonk_by_hand::verifier::Verifier;
use math::field;

fn main() {
    let mut srs = constants::SRS_BY_HAND;
    srs.generate_g_1_points();
    srs.generate_g_2_points();

    let field = constants::FIELD_17.clone();

    let pub_coin = constants::PUB_COIN.clone();

    const A: u32 = 3;
    const B: u32 = 4;
    const C: u32 = 5;

    let mut prover = Prover::new(field.clone(), vec![A, B, C], srs.copy());
    prover.set_public_coin(pub_coin.clone());
    let proof = prover.generate_proof();

    let mut verifier = Verifier::new(field, srs.copy());
    verifier.preprocess();
    verifier.provide_proof(pub_coin, proof);

    assert!(verifier.verify_proof());
}
