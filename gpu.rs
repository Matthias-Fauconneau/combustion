/*#![feature(format_args_capture)]#![allow(non_snake_case)]
mod vulkan;
fn wrap<const U: usize, const V: usize, const A: usize>(f: Function<U, V, A>) -> impl Fn([f64; U], [f64; V], [&[f64]; A]) -> Box<[f64]> { move |uniforms, values, arrays| {
	let values_arrays = Buffer::new(device, values.iter().copied())?;
	let output = Buffer::new(device, std::iter::repeat(0.).take(f.output))?;
	let buffers = [&[&load] as &[_], &[&store] as &[_]];
	pub fn as_bytes<T>(value: &T) -> &[u8] { unsafe{std::slice::from_raw_parts(value as *const T as *const u8, std::mem::size_of::<T>())} }
	let pipeline = device.pipeline(main.as_binary(), as_bytes(&uniforms), &buffers)?;
	device.bind(pipeline.descriptor_set, &buffers)?;
	let command_buffer = device.command_buffer(&pipeline, as_bytes(&uniforms), /*width*/1, states_len)?;
	let _time = device.submit_and_wait(command_buffer)?;
	store.map(device).unwrap()
}/*

fn main() -> Result<(), Box<dyn std::error::Error>> {
	/*let model = yaml_model::Loader::load_from_str(std::str::from_utf8(&std::fs::read(std::env::args().skip(1).next().unwrap())?)?)?;
	let model = yaml_model::parse(&model);
	use chemistry::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	let ref state = initial_state(&model);
	use {iter::map, itertools::Itertools};

	use vulkan::{Device, Buffer};
	let ref device = Device::new()?;

	if true {
		let rates = wrap(spirv::from(&reaction::rates(&species.thermodynamics[..species.len()-1], &map(&*model.reactions, |r| Reaction::new(species_names, r)))));
		assert!(state.volume == 1.);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let active_amounts = &amounts[0..amounts.len()-1];
		if let [energy_rate_RT, rates @ ..] = &*rates([*pressure_R], [total_amount, *T],[active_amounts]) {
			println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e3, (states_len as f32)/1e3/time);
			eprintln!("{}, HRR: {:.3e}", rates.iter().zip(&**species_names).format_with(", ", |(rate, name), f| f(&format!("{name}: {rate:.0}").to_string())), NA * kB * T * -energy_rate_RT);
		} else { unreachable!() }
	}
	#[cfg(feature="transport")] {
		let transport = wrap(transport::properties::<4>(&species));
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		if let [viscosity, thermal_conductivity, mixture_diffusion_coefficients @ ..] = &*transport([*pressure_R],[total_amount, *T],[amounts]) {
			eprintln!("μ: {viscosity:.4e}, λ: {thermal_conductivity:.4}, D: {:.4e}", mixture_diffusion_coefficients.iter().format(" "));
		} else { unreachable!() }
	}*/
	Ok(())
}
