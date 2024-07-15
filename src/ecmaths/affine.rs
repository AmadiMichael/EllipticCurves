use super::{jacobian::JacobianPoint, ru256::RU256};
use crate::curves::SECP256;
use primitive_types::U256;
use std::str::FromStr;

#[derive(Debug, Clone, PartialEq)]
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
