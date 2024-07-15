use crate::ecmaths::{affine::ECAffinePoint, ru256::RU256};

pub mod k1;
pub mod r1;

pub trait SECP256 {
    fn p() -> RU256;
    fn g() -> ECAffinePoint;
    fn n() -> RU256;
    fn a() -> RU256;
    fn b() -> RU256;
    fn n_div_2() -> RU256;
    fn sqrt_exp_num() -> RU256;
}
