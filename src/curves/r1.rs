use super::SECP256;
use crate::ecmaths::{affine::ECAffinePoint, ru256::RU256};
use primitive_types::U256;
use std::str::FromStr;

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

    fn sqrt_exp_num() -> RU256 {
        let p = Self::p();
        let p_plus_1 = p.v.checked_add(U256::one()).unwrap();
        return RU256 {
            v: p_plus_1.div_mod(U256::from(4)).0,
        };
    }
}
