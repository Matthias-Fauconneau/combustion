#![feature(thread_spawn_unchecked, unboxed_closures, fn_traits)]
use {fehler::throws, error::Error, combustion::{*, Property::*}};

#[throws] fn main() {
	let model = Model::new(model::Model::new(&std::fs::read("CH4+O2.ron")?)?);
	let (traps, function, rate) = model.rate::<{Volume}>();
	let handler = move |info:&libc::siginfo_t| {
		let offset = unsafe{(info.si_addr() as *const u8).offset_from(function as *const u8)};
		dbg!(traps.iter().find(|&x| dbg!(x.code_offset) == offset as u32).unwrap().source_location);
		std::process::abort()
	};
	unsafe{signal_hook_registry::register_unchecked(libc::SIGILL, handler)}?;
	let ref state = Simulation::new(&std::fs::read("CH4+O2.ron")?)?.state;
	let mut derivative = /*Derivative*/StateVector(std::iter::repeat(0.).take(model.len()).collect());
	println!("{:?}", {rate(state.constant(), &state.into(), &mut derivative); derivative});
}
