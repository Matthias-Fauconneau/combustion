//#![feature(default_free_fn, in_band_lifetimes, bindings_after_at, unboxed_closures)] #![allow(non_snake_case)]
/*
pub fn as_bytes<T>(value: &T) -> &[u8] { unsafe{std::slice::from_raw_parts(value as *const T as *const u8, std::mem::size_of::<T>())} }

fn time<T:Fn<()>>(task: T) -> (T::Output, f32) {
	let start = std::time::Instant::now();
	let result = task();
	(result, (std::time::Instant::now()-start).as_secs_f32())
}

fn print_time<T:Fn<()>>(task: T, header: impl std::fmt::Display) -> T::Output {
	let (result, time) = time(task);
	println!("{}: {:.1}ms", header, time*1e3);
	result
}

macro_rules! time { ($task:expr) => { print_time(|| { $task }, stringify!($task)) } }

mod vulkan;

use reaction::Simulation;*/
/*#[throws]*/ fn main() {
	let model = yaml_model::Loader::load_from_str(std::str::from_utf8(&std::fs::read(std::env::args().skip(1).next().unwrap())?)?)?;
	let model = yaml_model::parse(&model)?;
	use chemistry::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	let ref state = initial_state(&model);
	use {iter::map, itertools::Itertools, ast::wrap};
	if true {
		use reaction::*;
		let exp_Gibbs_RT = exp_Gibbs_RT(&species.thermodynamics[0..species.len()-1]);
	}
		let function = glsl(function)?;
	std::fs::write("/var/tmp/main.comp", &function)?;
  let main = shaderc::Compiler::new().unwrap().compile_into_spirv(&function, shaderc::ShaderKind::Compute, "main.comp", "main", None)?;

	use vulkan::{Device, Buffer};
	let ref device = Device::new()?;
	let load = Buffer::new(device, states.iter().copied())?;
	let store = Buffer::new(device, rates.iter().copied())?;
	let constants = reaction::Constants{
		pressure_Pa_R,
		temperature: 0,
		mass_fractions: (states_len*std::mem::size_of::<f64>()) as u32,
		mass_production_rates: (states_len*std::mem::size_of::<f64>()) as u32,
		heat_release_rate: 0,
		reference_temperature,
		mass_production_rate_factor,
		heat_release_rate_factor
	};
	let buffers = [&[&load] as &[_], &[&store] as &[_]];
	let pipeline = time!(device.pipeline(main.as_binary(), as_bytes(&constants), &buffers))?; // Compiles SPIRV -> Gen
	device.bind(pipeline.descriptor_set, &buffers)?;
	let width = 1;
	let command_buffer = device.command_buffer(&pipeline, as_bytes(&constants), width, states_len)?;
	dbg!();
	let time = device.submit_and_wait(command_buffer)?;
	dbg!();
	println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e3, (states_len as f32)/1e3/time);
	let rates = store.map(device).unwrap();
	reaction::report(&species_names, &rates);*/
}
