#[allow(dead_code)] pub const fn parse(s: &str) -> usize {
	let s = s.as_bytes();
	const POW10: [usize; 2] = [10, 1];
	let mut number: usize = 0;
	let mut index = 0;
	while index < s.len() {
		let digit = s[index] - b'0';
		//if digit > 9 { compile_error!("Invalid number"); }
		number += (digit as usize) * POW10[POW10.len()-s.len()+index];
		index += 1;
	}
	number
}
