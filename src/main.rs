use std::str::FromStr;

#[allow(unused_imports)]
use elliptic_curve_rust::{
    curves::{k1::K1, r1::R1, SECP256},
    ecmaths::ru256::RU256,
    signature::{PrivateKey, Signature},
};

fn main() {
    let curve = R1;
    let msg_hash = RU256::from_str("0x00").unwrap();
    let priv_key = PrivateKey::new(
        RU256::from_str("0xc1435991560e77992aaa190216c8939e3dc1855576a979963a3fd7110c04c316")
            .unwrap(),
    );

    let pub_key = priv_key.to_pub_key(&curve);
    println!("Public key: {:?}", &pub_key);

    let signature = priv_key.raw_sign(&msg_hash, &RU256::from_str("0x00").unwrap(), &curve);
    println!("Signature: {:?}", signature);

    println!(
        "Verify sig: {}",
        signature.raw_verify(&msg_hash, &pub_key, &curve)
    );

    println!(
        "Recover sig: {:?}",
        signature.raw_recover(&msg_hash, &curve)
    );
}
