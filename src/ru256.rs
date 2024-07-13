// modified from https://github.com/darth-cy/ecc/blob/main/src/ru256.rs
// temp solution for starting development

use crate::bytes;
use hex;
use primitive_types::U256;
use std::str::FromStr;

#[derive(Debug, Clone, PartialOrd)]
pub struct RU256 {
    pub v: U256,
}

#[derive(Debug, PartialEq, Eq)]
pub struct RU256ParseError;

impl FromStr for RU256 {
    type Err = RU256ParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match U256::from_str_radix(s, 16) {
            Ok(n) => return Ok(Self { v: n }),
            Err(_) => return Err(RU256ParseError),
        }
    }
}
impl ToString for RU256 {
    fn to_string(&self) -> String {
        let mut bytes = [0; 32];
        self.v.to_big_endian(&mut bytes);
        return hex::encode(bytes);
    }
}

impl RU256 {
    pub fn from_bytes(bs: &[u8]) -> Self {
        assert!(bs.len() <= 32, "big-endian");

        return Self {
            v: U256::from_big_endian(bs),
        };
    }
    pub fn to_bytes(&self, r: &mut [u8]) {
        self.v.to_big_endian(r);
    }

    // ******************************************************************
    // Identity Elements
    // 0 as additive identity; 1 as multiplicative identity;
    // ******************************************************************

    pub fn zero() -> Self {
        return Self::from_str("0x0").unwrap();
    }
    pub fn one() -> Self {
        return Self::from_str("0x1").unwrap();
    }

    // constant
    pub fn two() -> Self {
        return Self::from_str("0x2").unwrap();
    }
    pub fn three() -> Self {
        return Self::from_str("0x3").unwrap();
    }
    pub fn four() -> Self {
        return Self::from_str("0x4").unwrap();
    }
    pub fn eight() -> Self {
        return Self::from_str("0x8").unwrap();
    }

    // ******************************************************************
    // Modular Addition
    // Use: (a + b) mod p = ((a mod p) + (b mod p)) mod p;
    // ******************************************************************
    pub fn add_mod(&self, b: &Self, p: &Self) -> Self {
        let x1 = self.v.checked_rem(p.v).expect("modulo");
        let x2 = b.v.checked_rem(p.v).expect("modulo");

        let (mut x3, over) = x1.overflowing_add(x2);

        if over {
            x3 = x3
                .checked_add(
                    U256::MAX
                        .checked_sub(p.v)
                        .expect("sub")
                        .checked_add(U256::from_big_endian(&[1]))
                        .expect("convesion"),
                )
                .expect("add");
        }

        x3 = x3.checked_rem(p.v).expect("modulo");

        return Self { v: x3 };
    }

    // ******************************************************************
    // Modular Subtraction
    // Use: (a - b) mod p = ((a mod p) - (b mod p)) mod p;
    // Use: (-b) as congruent to (p - b) as additive inverse;
    // ******************************************************************
    pub fn sub_mod(&self, b: &Self, p: &Self) -> Self {
        let x1 = self.v.checked_rem(p.v).expect("modulo");
        let x2 = b.v.checked_rem(p.v).expect("modulo");

        return Self { v: x1 }.add_mod(&Self { v: (p.v - x2) }, p);
    }

    // ******************************************************************
    // Modular Multiplication
    // **NOTE: This is an add-and-double implementation.
    //         For performance, use primitive_types::U512 to hold product
    //         Cast to U256 after modulating in U512.
    //
    // Interpret multiplication as consecutive addition
    //
    // Use smaller multiplier (convert binary) as operational flags
    // ******************************************************************
    pub fn mul_mod(&self, b: &Self, p: &Self) -> Self {
        let x1 = Self {
            v: self.v.checked_rem(p.v).expect("modulo"),
        };
        let x2 = Self {
            v: b.v.checked_rem(p.v).expect("modulo"),
        };

        let mut base = Self::zero();

        let seq: Self;
        let adder: Self;

        if x1.v < x2.v {
            seq = x1.clone();
            adder = x2.clone();
        } else {
            seq = x2.clone();
            adder = x1.clone();
        }

        let mut seq_bytes = [0; 32];
        seq.to_bytes(&mut seq_bytes);

        let mut seq_binaries: Vec<u8> = vec![];
        bytes::bytes_to_binary(&seq_bytes, &mut seq_binaries);

        let mut on = false;
        for d in seq_binaries.into_iter() {
            if on {
                base = base.add_mod(&base, p);
            }
            if d > 0 {
                on = true;
                base = base.add_mod(&adder, p);
            }
        }

        return base;
    }

    // ******************************************************************
    // Modular Exponentiation
    // Interpret exponentiation as consecutive multiplication
    //
    // Use maller multiplier (convert binary) as operational flags
    // ******************************************************************
    pub fn exp_mod(&self, e: &Self, p: &Self) -> Self {
        let seq = e.clone();
        let multiplier = RU256 {
            v: self.v.checked_rem(p.v).expect("modulo"),
        };

        let mut base = Self::one();

        let mut seq_bytes = [0; 32];
        seq.to_bytes(&mut seq_bytes);

        let mut seq_binaries: Vec<u8> = vec![];
        bytes::bytes_to_binary(&seq_bytes, &mut seq_binaries);

        let mut on = false;
        for d in seq_binaries.into_iter() {
            if on {
                base = base.mul_mod(&base, p);
            }
            if d > 0 {
                on = true;
                base = base.mul_mod(&multiplier, p);
            }
        }

        return base;
    }

    // ******************************************************************
    // Modular Division
    // Use: (b^(-1)) as congruent to b^(p-2);
    //
    // Proof: b^(p-2) is the multiplicative inverse of b
    //        b^(p-2) * b = b^(p-1) = 1 (mod p)
    // ******************************************************************
    pub fn div_mod(&self, b: &Self, p: &Self) -> Self {
        assert!(p.v - 2 > U256::from_big_endian(&[0]));
        return self.mul_mod(&b.exp_mod(&RU256 { v: p.v - 2 }, p), p);
    }
}

impl PartialEq for RU256 {
    fn eq(&self, other: &Self) -> bool {
        return self.v == other.v;
    }
}

#[cfg(test)]
mod tests {
    use crate::ru256::RU256;
    use std::str::FromStr;

    #[test]
    fn ru256_addition_case_1() {
        let a = RU256::from_str("0xBD").unwrap();
        let b = RU256::from_str("0x2B").unwrap();
        let p = RU256::from_str("0xB").unwrap();

        let r = a.add_mod(&b, &p);

        assert_eq!(
            r.to_string(),
            "0000000000000000000000000000000000000000000000000000000000000001"
        );
    }

    #[test]
    fn ru256_addition_case_2() {
        let a = RU256::from_str("0xa167f055ff75c").unwrap();
        let b = RU256::from_str("0xacc457752e4ed").unwrap();
        let p = RU256::from_str("0xf9cd").unwrap();

        let r = a.add_mod(&b, &p);

        assert_eq!(
            r.to_string(),
            "0000000000000000000000000000000000000000000000000000000000006bb0"
        );
    }

    #[test]
    fn ru256_addition_case_3() {
        let a = RU256::from_str("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2E")
            .unwrap();
        let b = RU256::from_str("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2E")
            .unwrap();
        let p = RU256::from_str("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F")
            .unwrap();

        let r = a.add_mod(&b, &p);

        assert_eq!(
            r.to_string(),
            "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2d"
        );
    }

    #[test]
    fn ru256_subtraction_case_1() {
        let a = RU256::from_str("0x1ce606").unwrap(); // a = 189389.unwrap();
        let b = RU256::from_str("0xacc12484").unwrap(); // b = 289833894.unwrap();
        let p = RU256::from_str("0xf3fa3").unwrap(); // p = 99933.unwrap();

        let r = a.sub_mod(&b, &p);

        assert_eq!(
            r.to_string(),
            "000000000000000000000000000000000000000000000000000000000009645b"
        );
    }

    #[test]
    fn ru256_subtraction_case_2() {
        let a = RU256::from_str("0xacc12484").unwrap(); // a = 289833894.unwrap();
        let b = RU256::from_str("0x1ce606").unwrap(); // b = 189389.unwrap();
        let p = RU256::from_str("0xf3fa3").unwrap(); // p = 99933.unwrap();

        let r = a.sub_mod(&b, &p);

        assert_eq!(
            r.to_string(),
            "000000000000000000000000000000000000000000000000000000000005db48"
        );
    }

    #[test]
    fn ru256_multiplication_case() {
        let a = RU256::from_str("0xa167f055ff75c").unwrap(); // a = 283948457393954.unwrap();
        let b = RU256::from_str("0xacc457752e4ed").unwrap(); // b = 303934849383754.unwrap();
        let p = RU256::from_str("0xf9cd").unwrap(); // p = 6394.unwrap();

        let r = a.mul_mod(&b, &p);

        assert_eq!(
            r.to_string(),
            "000000000000000000000000000000000000000000000000000000000000e116"
        );
    }

    #[test]
    fn ru256_exponentiation_case() {
        let a = RU256::from_str("0x1ce606").unwrap(); // a = 189389.unwrap();
        let b = RU256::from_str("0xacc12484").unwrap(); // b = 289833894.unwrap();
        let p = RU256::from_str("0xf3fa3").unwrap(); // p = 99933.unwrap();

        let r = a.exp_mod(&b, &p);

        assert_eq!(
            r.to_string(),
            "000000000000000000000000000000000000000000000000000000000002a0fd"
        );
    }

    #[test]
    fn ru256_division_case() {
        let a = RU256::from_str("0x1ce606").unwrap(); // a = 189389.unwrap();
        let b = RU256::from_str("0xacc12484").unwrap(); // b = 289833894.unwrap();
        let p = RU256::from_str("0xf3fa3").unwrap(); // p = 99933.unwrap();

        let r = a.div_mod(&b, &p);

        assert_eq!(
            r.to_string(),
            "0000000000000000000000000000000000000000000000000000000000061f57"
        );
    }
}
