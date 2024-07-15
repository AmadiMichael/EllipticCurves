use super::{affine::ECAffinePoint, ru256::RU256};
use crate::curves::SECP256;
use primitive_types::U256;

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
