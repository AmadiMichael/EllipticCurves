use primitive_types::U256;

use crate::ru256::RU256;
use std::str::FromStr;

#[derive(Debug, Clone)]
pub struct ECAffinePoint {
    pub x: RU256,
    pub y: RU256,
}

impl ECAffinePoint {
    pub fn from_hex_coordinates(x: &str, y: &str) -> Self {
        return Self {
            x: RU256::from_str(x).unwrap(),
            y: RU256::from_str(y).unwrap(),
        };
    }
    pub fn to_hex_string(&self) -> String {
        return format!("04{}{}", self.x.to_string(), self.y.to_string());
    }
    pub fn is_zero_point(&self) -> bool {
        return self.x == RU256::zero() && self.y == RU256::zero();
    }

    // ******************************************************************
    // Identity Element
    // **NOTE: Imaginary. Implemented by setting both coordinates as 0
    //         need to check during operations
    // ******************************************************************

    pub fn zero_point() -> Self {
        return Self {
            x: RU256::zero(),
            y: RU256::zero(),
        };
    }

    pub fn add<T: SECP256>(&self, other: &Self, _: &T) -> Self {
        // checks
        assert!(self.y != other.y, "should use doubling");

        /*
         * Formula
         *
         *           y2 - y1
         * slope =  ---------
         *           x2 - x1
         *
         * x3 = slope ** 2 - x1 - x2
         * y3 = slope(x1 - x3) - y1
         */

        // implementation
        if self.is_zero_point() {
            return other.clone();
        }
        if other.is_zero_point() {
            return self.clone();
        }

        let p = &T::p();

        let slope = self
            .y
            .sub_mod(&other.y, p)
            .div_mod(&self.x.sub_mod(&other.x, p), p);

        let x = slope
            .mul_mod(&slope, p)
            .sub_mod(&self.x, p)
            .sub_mod(&other.x, p);
        let y = slope.mul_mod(&self.x.sub_mod(&x, p), p).sub_mod(&self.y, p);

        Self { x, y }
    }

    pub fn double<T: SECP256>(&self, _: &T) -> Self {
        /*
         * Formula
         *
         *           (3 * (x ** 2)) + a
         * slope =  -------------------
         *                2 * y
         *
         * x3 = slope ** 2 - x1 - x2
         * y3 = slope(x1 - x3) - y1
         */

        // implementation
        if self.is_zero_point() {
            return Self::zero_point();
        }
        if self.y == RU256::zero() {
            return Self::zero_point();
        }

        let p = &T::p();

        let slope = self
            .x
            .exp_mod(&RU256::two(), p)
            .mul_mod(&RU256::three(), p)
            .add_mod(&T::a(), p)
            .div_mod(&self.y.mul_mod(&RU256::two(), p), p);

        let x = slope
            .mul_mod(&slope, p)
            .sub_mod(&self.x, p)
            .sub_mod(&self.x, p);
        let y = slope.mul_mod(&self.x.sub_mod(&x, p), p).sub_mod(&self.y, p);

        Self { x, y }
    }

    pub fn multiply<T: SECP256>(&self, scalar: &RU256, curve: &T) -> Self {
        // Double and add method
        /*
         * R = 0
         * LOOP: R = (R * 2) + scalar_in_bit[i] * self
         * Note: i starts from 255 and goes down up until 0 (inclusive)
         */
        // implementation
        if self.y == RU256::zero() || scalar == &RU256::zero() {
            return Self::zero_point();
        }
        if scalar == &RU256::one() {
            return self.clone();
        }

        let mut r = Self::zero_point();

        let mut i = 255;
        while i != -1 {
            r = r.double(curve);
            let bit = (scalar.v >> i) & U256::one();
            if bit == U256::one() {
                r = r.add(&self, curve);
            }

            i -= 1;
        }

        r
    }

    pub fn to_jacobian(&self) -> JacobianPoint {
        JacobianPoint {
            x: self.x.clone(),
            y: self.y.clone(),
            z: RU256::one(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct JacobianPoint {
    pub x: RU256,
    pub y: RU256,
    pub z: RU256,
}

/**
 * Jacobian points add a third coordinate named z where
 * jacobian.x = affine.x / (z ** 2)
 * jacobian.y = affine.y / (z ** 3)
 */
impl JacobianPoint {
    pub fn is_zero_point(&self) -> bool {
        return self.x == RU256::zero() && self.y == RU256::zero();
    }

    pub fn zero_point() -> Self {
        return Self {
            x: RU256::zero(),
            y: RU256::zero(),
            z: RU256::one(),
        };
    }

    pub fn add<T: SECP256>(&self, other: &Self, curve: &T) -> Self {
        /*
         * u1 = x1 * (z2 ** 2)
         * u2 = x2 * (z1 ** 2)
         * s1 = y1 * (z2 ** 3)
         * s2 = y2 * (z1 ** 3)
         *
         * h = u2 - u1
         * r = s2 - s1
         *
         * h2 = h ** 2
         * h3 = h ** 3
         *
         * x3 = (r ** 2) - h3 - (2 * u1 * h2)
         * y3 = r * ((u1 * h2) - x3) - (s1 * h3)
         * z3 = h * z1 * z2
         */

        // implementation
        if self.is_zero_point() {
            return other.clone();
        }
        if other.is_zero_point() {
            return self.clone();
        }

        let p = &T::p();
        let z1z1 = self.z.mul_mod(&self.z, p);
        let z2z2 = other.z.mul_mod(&other.z, p);

        let u1 = self.x.mul_mod(&z2z2, p);
        let u2 = other.x.mul_mod(&z1z1, p);
        let s1 = self.y.mul_mod(&other.z.mul_mod(&z2z2, p), p);
        let s2 = other.y.mul_mod(&self.z.mul_mod(&z1z1, p), p);

        if u1 == u2 {
            if s1 != s2 {
                return Self::zero_point();
            }
            return self.double(curve);
        }

        let h = &u2.sub_mod(&u1, p);
        let h2 = &h.mul_mod(h, p);
        let h3 = &h2.mul_mod(h, p);

        let r = &s2.sub_mod(&s1, p);
        let v = &u1.mul_mod(h2, p);

        let x = r
            .mul_mod(&r, p)
            .sub_mod(h3, p)
            .sub_mod(&v.mul_mod(&RU256::two(), p), p);
        let y = r
            .mul_mod(&v.sub_mod(&x, p), p)
            .sub_mod(&s1.mul_mod(h3, p), p);
        let z = h.mul_mod(&self.z, p).mul_mod(&other.z, p);

        Self { x, y, z }
    }

    pub fn double<T: SECP256>(&self, _: &T) -> Self {
        /*
         * ysq = y ** 2
         * s = 4 * x * ysq
         * m = (3 * (x ** 2)) + (A * (z ** 4))
         *
         * x3 = (m ** 2) - 2 * s
         * y3 = m * (s - x3) - (8 * (ysq ** 2))
         * z3 = 2 * y1 * z1
         */

        if self.is_zero_point() {
            return Self::zero_point();
        }

        // implementation
        let p = &T::p();

        let ysq = self.y.mul_mod(&self.y, p);
        let s = self.x.mul_mod(&RU256::four(), p).mul_mod(&ysq, p);
        let m = self
            .x
            .mul_mod(&self.x, p)
            .mul_mod(&RU256::three(), p)
            .add_mod(&self.z.exp_mod(&RU256::four(), p).mul_mod(&T::a(), p), p);

        let x = m.mul_mod(&m, p).sub_mod(&RU256::two().mul_mod(&s, p), p);
        let y = m
            .mul_mod(&s.sub_mod(&x, p), p)
            .sub_mod(&RU256::eight().mul_mod(&ysq.mul_mod(&ysq, p), p), p);
        let z = self.y.mul_mod(&RU256::two(), p).mul_mod(&self.z, p);

        Self { x, y, z }
    }

    pub fn multiply<T: SECP256>(self, scalar: &RU256, curve: &T) -> Self {
        // Double and add method
        /*
         * R = 0
         * LOOP: R = (R * 2) + scalar_in_bit[i] * self
         * Note: i starts from 255 and goes down up until 0 (inclusive)
         */
        // implementation
        if self.y == RU256::zero() || scalar == &RU256::zero() {
            return Self::zero_point();
        }
        if scalar == &RU256::one() {
            return self;
        }

        let mut r = Self::zero_point();
        let mut i = 255;

        while i != -1 {
            r = r.double(curve);
            let bit = (scalar.v >> i) & U256::one();
            if bit == U256::one() {
                r = r.add(&self, curve);
            }

            i -= 1;
        }

        r
    }

    pub fn from_jacobian<T: SECP256>(&self, _: &T) -> ECAffinePoint {
        let p = &T::p();

        let z = RU256::one().div_mod(&self.z, p);
        let zz = z.mul_mod(&z, p);

        let x = self.x.mul_mod(&zz, p);
        let y = self.y.mul_mod(&zz.mul_mod(&z, p), p);

        ECAffinePoint { x, y }
    }
}

#[derive(Debug)]
pub struct Signature {
    pub r: RU256,
    pub s: RU256,
    pub v: RU256,
}

impl Signature {
    pub fn raw_sign<T: SECP256>(
        priv_key: &RU256,
        msg_hash: &RU256,
        nonce: &RU256,
        curve: &T,
    ) -> Self {
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
            .add_mod(&r.mul_mod(priv_key, n), n)
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

    // pub fn raw_recover(self, _msg_hash: RU256) -> ECAffinePoint {
    //     ECAffinePoint::zero_point()
    // }
}

pub trait SECP256 {
    fn p() -> RU256;
    fn g() -> ECAffinePoint;
    fn n() -> RU256;
    fn a() -> RU256;
    fn b() -> RU256;
    fn n_div_2() -> RU256;
}

#[derive(Debug)]
pub struct K1;

impl SECP256 for K1 {
    // ******************************************************************
    // SECP256K1 Curve Parameters
    // Reference: https://www.secg.org/sec2-v2.pdf
    // ******************************************************************

    fn p() -> RU256 {
        return RU256::from_str("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F")
            .unwrap();
    }
    fn g() -> ECAffinePoint {
        return ECAffinePoint {
            x: RU256::from_str("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798")
                .unwrap(),
            y: RU256::from_str("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8")
                .unwrap(),
        };
    }
    fn n() -> RU256 {
        return RU256::from_str("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141")
            .unwrap();
    }

    fn a() -> RU256 {
        RU256::zero()
    }

    fn b() -> RU256 {
        RU256::from_str("0x7").unwrap()
    }

    fn n_div_2() -> RU256 {
        RU256::from_str("0x7fffffffffffffffffffffffffffffff5d576e7357a4501ddfe92f46681b20a0")
            .unwrap()
    }
}

#[derive(Debug)]
pub struct R1;

impl SECP256 for R1 {
    // ******************************************************************
    // SECP256R1 Curve Parameters
    // Reference: https://www.secg.org/sec2-v2.pdf
    // ******************************************************************

    fn p() -> RU256 {
        return RU256::from_str(
            "0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF",
        )
        .unwrap();
    }
    fn g() -> ECAffinePoint {
        return ECAffinePoint {
            x: RU256::from_str(
                "0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296",
            )
            .unwrap(),
            y: RU256::from_str(
                "0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5",
            )
            .unwrap(),
        };
    }
    fn n() -> RU256 {
        return RU256::from_str(
            "0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551",
        )
        .unwrap();
    }

    fn a() -> RU256 {
        RU256::from_str("0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC")
            .unwrap()
    }

    fn b() -> RU256 {
        RU256::from_str("0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B")
            .unwrap()
    }

    fn n_div_2() -> RU256 {
        RU256::from_str("0x7fffffff800000007fffffffffffffffde737d56d38bcf4279dce5617e3192a8")
            .unwrap()
    }
}
