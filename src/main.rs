#![allow(mixed_script_confusables, non_snake_case, incomplete_features, confusable_idents)]
#![feature(type_ascription, array_map, non_ascii_idents, const_generics, const_evaluatable_checked, destructuring_assignment, test)]
use {fehler::throws, error::Error, combustion::{*, Property::*}};

#[throws] fn main() {
	pretty_env_logger::init();
	let model = Model::new(model::Model::new(&std::fs::read("CH4+O2.ron")?)?);
	let rate = model.rate::<{Volume}>();
	let ref state = Simulation::new(&std::fs::read("CH4+O2.ron")?)?.state;
	let mut derivative = /*Derivative*/StateVector(std::iter::repeat(0.).take(model.len()).collect());
	println!("{:?}", {rate(state.constant(), &state.into(), &mut derivative); &derivative});
	/*let len = 100000;
	let constant = state.constant();
	let state = state.into();
	let start = std::time::Instant::now();
	for _ in 0..len { rate(constant, &state, &mut derivative); }
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f32();
	println!("{:.1}ms\t{:.0}K/s", time*1e3, (len as f32)/time/1e3);*/
}
