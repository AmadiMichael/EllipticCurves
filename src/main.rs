use std::str::FromStr;

use crate::{
    ru256::RU256,
    secp256::{JacobianPoint, Signature, K1, R1, SECP256},
};

pub mod bytes;
pub mod ru256;
pub mod secp256;

fn main() {
    // let point1 = secp256::ECAffinePoint::from_hex_coordinates(
    //     "0x7cf27b188d034f7e8a52380304b51ac3c08969e277f21b35a60b48fc47669978",
    //     "0x07775510db8ed040293d9ac69f7430dbba7dade63ce982299e04b79d227873d1",
    // );
    // let point2 = secp256::ECAffinePoint::from_hex_coordinates(
    //     "0x5ecbe4d1a6330a44c8f7ef951d4bf165e6c6b721efada985fb41661bc6e7fd6c",
    //     "0x8734640c4998ff7e374b06ce1a64a2ecd82ab036384fb83d9a79b127a27d5032",
    // );

    // println!("{:?}", point1.add(&point2, R1));
    // println!(
    //     "{:?}",
    //     JacobianPoint::from_jacobian(
    //         &JacobianPoint::add(
    //             &JacobianPoint::to_jacobian(point1),
    //             &JacobianPoint::to_jacobian(point2),
    //             &R1
    //         ),
    //         &R1
    //     )
    // );

    // let point1 = secp256::ECAffinePoint::from_hex_coordinates(
    //     "0x7cf27b188d034f7e8a52380304b51ac3c08969e277f21b35a60b48fc47669978",
    //     "0x07775510db8ed040293d9ac69f7430dbba7dade63ce982299e04b79d227873d1",
    // );
    // println!("{:?}", point1.double(R1));
    // println!(
    //     "{:?}",
    //     JacobianPoint::from_jacobian(&JacobianPoint::to_jacobian(point1).double(&R1), &R1)
    // );

    println!(
        "Public key: {:?}",
        JacobianPoint::to_jacobian(K1::g())
            .multiply(
                &RU256::from_str(
                    "0xc2362fa72dcc0e1201473ee43071b1750e4879bb67bdd75e7d7efe963aa2d5a9"
                )
                .unwrap(),
                &K1
            )
            .from_jacobian(&K1)
    );

    println!(
        "Signature: {:?}",
        Signature::raw_sign(
            &RU256::from_str("0xc2362fa72dcc0e1201473ee43071b1750e4879bb67bdd75e7d7efe963aa2d5a9")
                .unwrap(),
            &RU256::from_str("0x04").unwrap(),
            &RU256::from_str("0x03").unwrap(),
            &K1
        )
    )
}
