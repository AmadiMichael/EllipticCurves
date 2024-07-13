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

    pub fn zero_point() -> ECAffinePoint {
        return ECAffinePoint {
            x: RU256::zero(),
            y: RU256::zero(),
        };
    }

    pub fn add<T: SECP256>(&self, other: &Self, _: T) -> Self {
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

    pub fn double<T: SECP256>(&self, _: T) -> Self {
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
            if s1 == s2 {
                return Self {
                    x: RU256::zero(),
                    y: RU256::zero(),
                    z: RU256::two(),
                };
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
         * m = (3 * (x ** 2)) + (A * z1)
         *
         * x3 = (m ** 2) - 2 * s
         * y3 = m * (s - x3) - (8 * (ysq ** 2))
         * z3 = 2 * y1 * z1
         */

        // implementation
        let p = &T::p();

        let ysq = self.y.mul_mod(&self.y, p);
        let s = self.x.mul_mod(&RU256::four(), p).mul_mod(&ysq, p);
        let m = self
            .x
            .mul_mod(&self.x, p)
            .mul_mod(&RU256::three(), p)
            .add_mod(&self.z.mul_mod(&T::a(), p), p);

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
         * Note: i starts from 256 and goes down
         */
        // implementation
        if self.y == RU256::zero() || scalar == &RU256::zero() {
            return Self {
                x: RU256::zero(),
                y: RU256::zero(),
                z: RU256::one(),
            };
        }
        if scalar == &RU256::one() {
            return self;
        }

        let mut r = JacobianPoint::to_jacobian(T::zero());

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

    pub fn to_jacobian(aff: ECAffinePoint) -> Self {
        Self {
            x: aff.x,
            y: aff.y,
            z: RU256::one(),
        }
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
    ) -> Signature {
        let encoded_nonce = JacobianPoint::from_jacobian(
            &JacobianPoint::to_jacobian(T::g()).multiply(nonce, curve),
            curve,
        );

        let n = &T::n();

        let r = encoded_nonce.x;
        let mut s = msg_hash
            .add_mod(&r.mul_mod(priv_key, n), n)
            .div_mod(&nonce, n);
        let mut v = RU256::from_str("0x1b").unwrap();

        // use lower order of n
        if s < n.div_mod(&RU256::two(), n) {
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

    pub fn raw_recover(self, _msg_hash: RU256) -> ECAffinePoint {
        ECAffinePoint::zero_point()
    }
}

pub trait SECP256 {
    fn p() -> RU256;
    fn g() -> ECAffinePoint;
    fn n() -> RU256;
    fn a() -> RU256;
    fn b() -> RU256;
    fn zero() -> ECAffinePoint;
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

    fn zero() -> ECAffinePoint {
        ECAffinePoint {
            x: RU256::zero(),
            y: RU256::zero(),
        }
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
        return RU256::from_str("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff")
            .unwrap();
    }
    fn g() -> ECAffinePoint {
        return ECAffinePoint {
            x: RU256::from_str("6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296")
                .unwrap(),
            y: RU256::from_str("4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5")
                .unwrap(),
        };
    }
    fn n() -> RU256 {
        return RU256::from_str("ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551")
            .unwrap();
    }

    fn a() -> RU256 {
        RU256::from_str("ffffffff00000001000000000000000000000000fffffffffffffffffffffffc").unwrap()
    }

    fn b() -> RU256 {
        RU256::from_str("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b").unwrap()
    }

    fn zero() -> ECAffinePoint {
        ECAffinePoint {
            x: RU256::zero(),
            y: RU256::zero(),
        }
    }
}
