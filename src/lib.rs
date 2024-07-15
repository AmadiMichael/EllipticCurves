pub mod bytes;
pub mod curves;
pub mod ecmaths;
pub mod signature;

#[cfg(test)]
mod tests {
    use crate::{curves, ecmaths::ru256::RU256, signature::PrivateKey};
    use std::str::FromStr;

    #[test]
    fn test_ec_sign_and_verify() {
        let curve = curves::k1::K1;

        let priv_key = PrivateKey::new(
            RU256::from_str("0xc1435991560e77992aaa190216c8939e3dc1855576a979963a3fd7110c04c316")
                .unwrap(),
        );
        let pub_key = priv_key.to_pub_key(&curve);

        for i in 0..10 {
            let msg_hash = RU256::from_str(format!("0x0{}", i).as_str()).unwrap();
            let nonce = RU256::from_str(format!("0x0{}", 10 - i).as_str()).unwrap();

            let signature = priv_key.raw_sign(&msg_hash, &nonce, &curve);

            assert!(
                signature.raw_verify(&msg_hash, &pub_key, &curve),
                "raw verify failed in iteration {}",
                i
            );
            assert_eq!(
                signature.raw_recover(&msg_hash, &curve),
                pub_key,
                "raw recover failed in iteration {}",
                i
            );
        }
    }
}
