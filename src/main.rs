#![feature(type_ascription)]#![feature(non_ascii_idents)]#![allow(confusable_idents,non_snake_case,unused_variables,unused_mut)]
use {fehler::throws, error::Error, combustion::{*, Property::*}};

#[throws] fn main() {
	let model = &std::fs::read("CH4+O2.ron")?;
	let model = model::Model::new(&model)?;
	let ref state = Simulation::new(&model)?.state;
	let model = Model::new(model);
	let (_traps, (_function, _size), rate) = model.rate::<{Volume}>();
	let mut derivative = /*Derivative*/StateVector::<{Volume}>(std::iter::repeat(0.).take(/*model.len()*/2+model.reactions.len()).collect());
	/*{
		let function = unsafe{std::slice::from_raw_parts(function as *const u8, size)}; // refcheck: leaked from dropped JITModule
		let (constant, ref state, derivative) = (state.constant::<{Volume}>(), state.into(): StateVector::<{Volume}>, &mut derivative);
		let constant = constant.0 as f32;
		use x86emu::*;
		let mut guest = State::new();
		allocate_stack(&mut guest);
		load(&mut guest, function);
		let mut heap = Heap::new(&mut guest);
		let state = heap.push_slice(&mut guest, &state.0);
		let derivative = heap.push_slice(&mut guest, &derivative.0);
		call(&mut guest, &[state as i64, derivative as i64], &[constant]);
		pretty_env_logger::init();
		guest.print_instructions = true;
		guest.execute()
	}*/
	let rate = {rate(state.constant(), &state.into(), &mut derivative); derivative.0};
	for (i, v) in rate.iter().enumerate() { if v.abs() > 1e-29 { print!("{}:{:.3e} ", i, v); } }
	println!("");
}
