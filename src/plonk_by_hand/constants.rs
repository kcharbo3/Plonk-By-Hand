use crate::{CurvePoint, ExtensionCurvePoint, Field, PublicCoin, ECC, SRS};

pub const FIELD_17: Field = Field { order: 17 };

pub const SRS_BY_HAND: SRS = SRS {
    g_1: CurvePoint { x: 1, y: 2 },
    g_2: ExtensionCurvePoint {
        x: 36,
        y: 31,
        u: true,
    },
    g_1_points: Vec::new(),
    g_2_points: Vec::new(),
    degree: 4,
    s: 2,
    scalar_field: FIELD_17,
    ecc: ECC {
        field: Field { order: 101 },
    },
};

pub const PUB_COIN: PublicCoin = PublicCoin {
    b1: 7,
    b2: 4,
    b3: 11,
    b4: 12,
    b5: 16,
    b6: 2,
    b7: 14,
    b8: 11,
    b9: 7,
    alpha: 15,
    beta: 12,
    gamma: 13,
    zed: 5,
    v: 12,
    u: 4,
};
