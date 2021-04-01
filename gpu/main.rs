#![feature(default_free_fn, in_band_lifetimes, bindings_after_at, unboxed_closures)] #![allow(non_snake_case)]
mod vulkan;

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

fn benchmark<T>(task: impl Fn()->T, times: usize, header: impl std::fmt::Display) -> T {
	let result = task();
	println!("{}: {:.1}ms", header, time(|| for _ in 0..times { task(); }).1/(times as f32)*1e3);
	result
}

macro_rules! benchmark { ($task:expr, $times:expr) => { benchmark(|| { $task }, $times, stringify!($task)) } }

#[fehler::throws(Box<dyn std::error::Error>)] fn main() {
	use iter::box_collect;
	let model = &std::fs::read("CH4+O2.ron")?;
	use combustion::{*, transport::*};
	let model = model::Model::new(&model)?;
	let ref state = Simulation::new(&model)?.state;
	//#[cfg(feature="transport")] {
	let (_species_names, species) = Species::new(model.species);
	let ref transport_polynomials = species.transport_polynomials();
	let transport = benchmark!(transport::transport(&species.molar_mass, transport_polynomials, state), 1);

	use vulkan::{Device, Buffer};
	let ref device = Device::new()?;
	let stride = 32;
	let len = 32/stride*stride;
	let State{temperature, pressure, amounts, ..} = state;
	let temperature = Buffer::new(device, (0..len).map(|_| temperature))?;
	let amounts_buffers = box_collect(amounts.iter().map(|n| Buffer::new(device, (0..len).map(|_| n)).unwrap()));
	let ref_amounts = box_collect(amounts_buffers.iter().map(|n| n));
	let d_temperature = Buffer::new(device, (0..len).map(|_| f64::NAN))?;
	let d_amounts = box_collect(amounts.iter().map(|_| Buffer::new(device, (0..len).map(|_| f64::NAN)).unwrap()));
	let ref_d_amounts = box_collect(d_amounts.iter().map(|n| n));
	let viscosity = Buffer::new(device, (0..len).map(|_| f64::NAN))?;
	let thermal_conductivity = Buffer::new(device, (0..len).map(|_| f64::NAN))?;
	let mixture_averaged_thermal_diffusion_coefficients = iter::eval(len, |_| Buffer::new(device, (0..len).map(|_| f64::NAN)).unwrap());
	let ref_mixture_averaged_thermal_diffusion_coefficients = box_collect(mixture_averaged_thermal_diffusion_coefficients.iter().map(|n| n));
	let ref buffers = [&[&temperature] as &[_], &ref_amounts, &[&d_temperature] as &[_], &ref_d_amounts, &[&viscosity] as &[_], &[&thermal_conductivity] as &[_], &ref_mixture_averaged_thermal_diffusion_coefficients];
	let ref constants = [pressure];
	let ref pipeline = time!(device.pipeline(constants, buffers))?; // Compiles SPIRV -> Gen
	device.bind(pipeline.descriptor_set, buffers)?;
	let command_buffer  = device.command_buffer(pipeline, constants, stride, len)?;
	for _ in 0..1 {
		let time = device.submit_and_wait(command_buffer)?;
		let all_same = |buffer:&Buffer| {
			let buffer = buffer.map(device).unwrap();
			assert!(buffer.len() == len);
			for &v in buffer.iter() { assert_eq!(v, buffer[0]); }
			buffer[0]
		};
		let gpu_transport = Transport{
			viscosity: all_same(&viscosity),
			thermal_conductivity: all_same(&thermal_conductivity),
			mixture_averaged_thermal_diffusion_coefficients: box_collect(mixture_averaged_thermal_diffusion_coefficients.iter().map(|buffer| all_same(buffer))),
		};
		use AbsError;
		if transport.error(&gpu_transport) > 3e-6 { println!("{:?}\n{:?}", transport, gpu_transport); }
		println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", len as f32/1e3, time*1e3, time/(len as f32)*1e3, (len as f32)/1e3/time);
	}
}
