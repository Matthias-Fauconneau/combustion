#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds)]#![allow(non_snake_case)]
mod vulkan;

use {iter::{box_, map}, ast::*, vulkan::*, fehler::throws, anyhow::{Error,Result}};
#[throws] fn wrap<const U: usize>(device: &'t Device, function: &Function) -> impl 't+Fn([f64; U], &[&[f64]]) -> Result<Box<[Box<[f64]>]>> {
	let input_len = function.input;
	let output_len = function.output;
	let function = spirv::from(U, function)?;
	move |uniforms:[f64; U], input:&[&[f64]]| {
		assert!(U+input.len() == input_len);
		let states_len = input[0].len();
		let input = map(&*input, |array| Buffer::new(device, array.iter().copied()).unwrap());
		let output = map(0..output_len, |_| Buffer::new(device, vec![0.; states_len]).unwrap());
		let buffers = box_(input.iter().chain(&*output));
		pub fn as_bytes<T>(value: &T) -> &[u8] { unsafe{std::slice::from_raw_parts(value as *const T as *const u8, std::mem::size_of::<T>())} }
		let pipeline = device.pipeline(&function, as_bytes(&uniforms), &buffers)?;
		device.bind(pipeline.descriptor_set, &buffers)?;
		let command_buffer = device.command_buffer(&pipeline, as_bytes(&uniforms), /*width*/1, states_len)?;
		let _time = device.submit_and_wait(command_buffer)?;
		Ok(map(&*output, |array| (*array.map(device).unwrap()).into()))
	}
}

#[throws] fn main() {
	let model = yaml_model::Loader::load_from_str(std::str::from_utf8(&std::fs::read(std::env::args().skip(1).next().unwrap())?)?)?;
	let model = yaml_model::parse(&model);
	use chemistry::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	let ref state = initial_state(&model);
	use itertools::Itertools;
	let ref device = Device::new()?;
	if true {
		let rates = wrap(device, &reaction::rates(&species.thermodynamics[..species.len()-1], &map(&*model.reactions, |r| Reaction::new(species_names, r))))?;
		assert!(state.volume == 1.);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let states_len = 1;
		let active_amounts = map(&amounts[0..amounts.len()-1], |&n| vec![n; states_len].into_boxed_slice());
		let_!{ [energy_rate_RT, rates @ ..] = &*rates([*pressure_R], &*([&[&[total_amount] as &[_], &[*T]] as &[_], &*map(&*active_amounts, |a| &**a)].concat()))? => {
		//println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e3, (states_len as f32)/1e3/time);
		let energy_rate_RT = energy_rate_RT[0];
		let rates = rates.iter().map(|a| a[0]);
		eprintln!("{}, HRR: {:.3e}", rates.zip(&**species_names).format_with(", ", |(rate, name), f| f(&format!("{name}: {rate:.0}").to_string())), NA * kB * T * -energy_rate_RT);
		}}
	}
	#[cfg(feature="transport")] {
		let transport = wrap(device, transport::properties::<4>(&species));
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let_!{ [viscosity, thermal_conductivity, mixture_diffusion_coefficients @ ..] = &*transport([*pressure_R], &[&[total_amount, *T], active_amounts].concat()) => {
			eprintln!("μ: {viscosity:.4e}, λ: {thermal_conductivity:.4}, D: {:.4e}", mixture_diffusion_coefficients.iter().format(" "));
		}}
	}
}