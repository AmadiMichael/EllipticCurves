use super::SECP256;
use crate::ecmaths::{affine::ECAffinePoint, ru256::RU256};
use primitive_types::U256;
use std::str::FromStr;

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

    fn sqrt_exp_num() -> RU256 {
        let p = Self::p();
        let p_plus_1 = p.v.checked_add(U256::one()).unwrap();
        return RU256 {
            v: p_plus_1.div_mod(U256::from(4)).0,
        };
    }
}
