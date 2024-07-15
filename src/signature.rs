use crate::{curves::SECP256, ecmaths::affine::ECAffinePoint, ecmaths::ru256::RU256};
use primitive_types::U256;
use std::str::FromStr;

pub struct PrivateKey(RU256);
impl PrivateKey {
    pub fn new(key: RU256) -> Self {
        Self(key)
    }

    pub fn to_pub_key<T: SECP256>(&self, curve: &T) -> ECAffinePoint {
        return T::g()
            .to_jacobian()
            .multiply(&self.0, curve)
            .from_jacobian(curve);
    }

    pub fn raw_sign<T: SECP256>(&self, msg_hash: &RU256, nonce: &RU256, curve: &T) -> Signature {
        /*
         * k = nonce
         * r, y = (k * G).x, (k * G).y
         * s = 1/k * (h + (p * r))
         * v = 27 + xor((s < half_n), (y % 2 == 0))
         */
        let n = &T::n();

        let encoded_nonce = T::g()
            .to_jacobian()
            .multiply(nonce, curve)
            .from_jacobian(curve);
        let r = encoded_nonce.x;
        let mut s = msg_hash
            .add_mod(&r.mul_mod(&self.0, n), n)
            .div_mod(&nonce, n);
        let mut v = RU256::from_str("0x1b").unwrap();

        // use lower order of n
        if s <= T::n_div_2() {
            v = match encoded_nonce.y.v % 2 == U256::zero() {
                true => v,
                false => v.add_mod(&RU256::one(), n),
            }
        } else {
            s = n.sub_mod(&s, n);
            v = match encoded_nonce.y.v % 2 == U256::zero() {
                true => v.add_mod(&RU256::one(), n),
                false => v,
            }
        }

        Signature { r, s, v }
    }
}

#[derive(Debug)]
pub struct Signature {
    pub r: RU256,
    pub s: RU256,
    pub v: RU256,
}

impl Signature {
    pub fn raw_verify<T: SECP256>(
        &self,
        msg_hash: &RU256,
        pub_key: &ECAffinePoint,
        curve: &T,
    ) -> bool {
        /*
         * sInv = 1/s
         * a = G * (sInv * h)
         * b = PubKey * (sInv * r)
         * c = a + b
         * c.x == r
         */

        let n = &T::n();

        let a = &T::g()
            .to_jacobian()
            .multiply(&msg_hash.div_mod(&self.s, n), curve);
        let b = pub_key
            .to_jacobian()
            .multiply(&self.r.div_mod(&self.s, n), curve);
        let c = a.add(&b, curve);

        return c.from_jacobian(curve).x == self.r;
    }

    pub fn raw_recover<T: SECP256>(self, _msg_hash: &RU256, curve: &T) -> ECAffinePoint {
        /*
         * assert that x is a valid point on curve
         *
         * PubKey = ((encoded_nonce * s) + (G * (-h))) / r
         */

        assert!(
            self.v == RU256::from_str("0x1b").unwrap()
                || self.v == RU256::from_str("0x1c").unwrap(),
            "invalid V",
        );

        let p = &T::p();
        let n = &T::n();

        // prove that self.r is a valid x on elliptic curve y**2 = x**3 + ax + b
        let x_cubed_ax_b = self
            .r
            .exp_mod(&RU256::three(), p)
            .add_mod(&T::a().mul_mod(&self.r, p), p)
            .add_mod(&T::b(), p);
        let possible_y = x_cubed_ax_b.exp_mod(&T::sqrt_exp_num(), p);
        let y = match (self.v.v.div_mod(U256::from(2)).1 == U256::one())
            ^ (possible_y.v.div_mod(U256::from(2)).1 == U256::one())
        {
            true => possible_y,
            false => T::p().sub_mod(&possible_y, p),
        };

        assert_eq!(
            x_cubed_ax_b.sub_mod(&y.mul_mod(&y, p), p),
            RU256::zero(),
            "sig invalid, r cannot be x coordinate of a point of the curve",
        );
        assert!(
            self.r.v.div_mod(T::n().v).1 != U256::zero()
                && self.s.v.div_mod(T::n().v).1 != U256::zero(),
            "r % n or s % n is 0"
        );

        let a = ECAffinePoint {
            x: self.r.clone(),
            y,
        }
        .to_jacobian()
        .multiply(&self.s, curve);
        let b = T::g()
            .to_jacobian()
            .multiply(&n.sub_mod(&_msg_hash, n), curve);
        let c = a.add(&b, curve);
        let pub_key = c.multiply(&RU256::one().div_mod(&self.r, n), curve);

        pub_key.from_jacobian(curve)
    }
}
