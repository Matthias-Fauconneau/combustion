#![feature(type_ascription)]#![feature(non_ascii_idents)]#![allow(confusable_idents,non_snake_case,unused_variables,unused_mut)]
use {fehler::throws, error::Error, combustion::*};

#[throws] fn main() {
	let model = &std::fs::read("CH4+O2.ron")?;
	let model = model::Model::new(&model)?;
	let ref state = Simulation::new(&model)?.state;
	#[cfg(feature="transport")] {
	}
	#[cfg(feature="reaction")] {
		let (_, rate) = model.rate();
		let mut derivative = /*Derivative*/StateVector::<{Property::Volume}>(std::iter::repeat(0.).take(2+model.len()-1).collect());

		{
			let rate = {rate(state.constant(), &state.into(), &mut derivative); &derivative.0};
			for (i, v) in rate.iter().enumerate() { if v.abs() > 1e-29 { print!("{}:{:.3e} ", i, v); } }
			println!("");
		}
		let len = 100000;
		let constant = state.constant();
		let state = state.into();
		let start = std::time::Instant::now();
		for _ in 0..len { rate(constant, &state, &mut derivative) }
		let end = std::time::Instant::now();
		let time = (end-start).as_secs_f32();
		println!("{:.1}ms\t{:.0}K/s", time*1e3, (len as f32)/time/1e3);
	}
}
