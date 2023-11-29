pub fn bytes_to_binary(i: &[u8; 32], r: &mut Vec<u8>) {
    for m in i.iter() {
        format!("{:8b}", m).chars().for_each(|b| {
            if b == '1' {
                r.push(1);
            } else {
                r.push(0)
            }
        });
    }
}
