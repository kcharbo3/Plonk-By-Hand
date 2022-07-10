// Note: Does not handle overflows
#[derive(Clone)]
pub struct Field {
    pub order: u32,
}

impl Field {
    pub fn add(&self, a: u32, b: u32) -> u32 {
        (a + b) % self.order
    }

    pub fn subtract(&self, a: u32, b: u32) -> u32 {
        self.add(a, self.additive_inverse(b))
    }

    pub fn multiply(&self, a: u32, b: u32) -> u32 {
        (a * b) % self.order
    }

    pub fn divide(&self, a: u32, b: u32) -> u32 {
        self.multiply(a, self.multiplicative_inverse(b))
    }

    pub fn exponent(&self, a: u32, b: u32) -> u32 {
        if b == 0 {
            return 1;
        }

        let mut res = a;
        for _i in 1..b {
            res = self.multiply(res, a);
        }

        res % self.order
    }

    pub fn additive_inverse(&self, a: u32) -> u32 {
        (self.order - (a % self.order)) % self.order
    }

    pub fn multiplicative_inverse(&self, a: u32) -> u32 {
        // TODO: change to euclidean algo
        let mut inv = 0;
        for i in 0..self.order {
            if self.multiply(i, a) == 1 {
                inv = i;
                break;
            }
        }
        inv
    }
}
