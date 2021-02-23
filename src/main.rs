#![allow(mixed_script_confusables, non_snake_case, incomplete_features, confusable_idents)]
#![feature(type_ascription, array_map, non_ascii_idents, const_generics, const_evaluatable_checked, destructuring_assignment, test)]
use {fehler::throws, ron::Error, combustion::{*, Property::*}};

#[throws] fn main() {
	let system = convert::system!("CH4+O2.ron");
	let state = Simulation::new(include_str!("../CH4+O2.ron"))?.state;
	println!("{:?}", system.rate::<{Volume}>(&state));
	let len = 10000;
	let start = std::time::Instant::now();
	for _ in 0..len { system.rate::<{Volume}>(&state); }
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f32();
	println!("{:.1}ms\t{:.0}K/s", time*1e3, (len as f32)/time/1e3);
}
