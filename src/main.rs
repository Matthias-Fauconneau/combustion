#![allow(mixed_script_confusables, non_snake_case, incomplete_features, confusable_idents)]
#![feature(type_ascription, array_map, non_ascii_idents, const_generics, const_evaluatable_checked, destructuring_assignment, test)]
use {fehler::throws, ron::Error, model::Troe, system::{NASA7, RateConstant, Property::{self, *}}, combustion::Simulation};

r#macro::r#macro!{rate_and_jacobian}

#[throws] fn main() {
	let ref state = Simulation::new(include_str!(concat!("../",concat!(env!("MODEL"),".ron"))))?.state;
	const CONSTANT: Property = Volume;
	println!("{:?}", rate_and_jacobian::<CONSTANT>(state.constant::<CONSTANT>(), &state.into()).unwrap().0);
	let len = 100000;
	let constant = state.constant::<CONSTANT>();
	let ref state = state.into();
	let start = std::time::Instant::now();
	for _ in 0..len { rate_and_jacobian::<CONSTANT>(constant, state).unwrap().0; }
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f32();
	println!("{:.1}ms\t{:.0}K/s", time*1e3, (len as f32)/time/1e3);
}
