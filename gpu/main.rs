#![feature(default_free_fn, in_band_lifetimes, bindings_after_at)] #![allow(non_snake_case)]
mod vulkan;

#[fehler::throws(anyhow::Error)] fn main() {
	use iter::{vec::{eval, generate, Vector}, box_collect};
	let system = std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let Simulation{system, pressure_R, state: state@combustion::State{temperature, amounts}, ..} = Simulation::<35>::new(&system)?;
	//let f = system.dt_J(pressure_r, &{use iter::into::IntoChain; from_iter([temperature].chain(amounts))});
	//let f = (f.0, /*f.1.concat().try_into().unwrap()*/);
	let transport = system.transport(pressure_R, &state);

	use vulkan::{Device, Buffer};
	let ref device = Device::new()?;
	let stride = 1;
	let len = 1/stride*stride;
	let temperature = Buffer::new(device, (0..len).map(|_| temperature))?;
	let amounts_buffers = eval(amounts, |n| Buffer::new(device, (0..len).map(|_| n)).unwrap());
	let ref_amounts = box_collect(amounts_buffers.iter().map(|n| n));
	let d_temperature = Buffer::new(device, (0..len).map(|_| f64::NAN))?;
	let d_amounts = eval(amounts, |_| Buffer::new(device, (0..len).map(|_| f64::NAN)).unwrap());
	let ref_d_amounts = box_collect(d_amounts.iter().map(|n| n));
	let viscosity = Buffer::new(device, (0..len).map(|_| f64::NAN))?;
	let thermal_conductivity = Buffer::new(device, (0..len).map(|_| f64::NAN))?;
	let mixture_averaged_thermal_diffusion_coefficients = generate(|_| Buffer::new(device, (0..len).map(|_| f64::NAN)).unwrap()).collect();
	let ref_mixture_averaged_thermal_diffusion_coefficients = box_collect(mixture_averaged_thermal_diffusion_coefficients.iter().map(|n| n));
	let ref buffers = [&[&temperature] as &[_], &ref_amounts, &[&d_temperature] as &[_], &ref_d_amounts/*, &box_collect(jacobian.iter_mut().map(|n| n))*/, &[&viscosity] as &[_], &[&thermal_conductivity] as &[_], &ref_mixture_averaged_thermal_diffusion_coefficients];
	//let mut jacobian: [_;(2+S-1)*(2+S-1)] = generate(|_| Buffer::new(device, (0..len).map(|_| f64::NAN)).unwrap());
	let ref constants = [pressure_R];
	let ref pipeline = device.pipeline(constants, buffers)?;
	device.bind(pipeline.descriptor_set, buffers)?;
	let command_buffer  = device.command_buffer(pipeline, constants, stride, len)?;
	//let _warmup = device.submit_and_wait(command_buffer)?;
	for _ in 0..1 {
		let start = std::time::Instant::now();
		let time = device.submit_and_wait(command_buffer)?;
		let end = std::time::Instant::now();
		let cpu_time = (end-start).as_secs_f32();
		assert!(cpu_time > time);
		/*let gpu_f = (
			*(box_collect(std::iter::once(&d_temperature).chain(d_amounts.iter()).map(|buffer:&Buffer| {
				let buffer = buffer.map(device).unwrap();
				assert!(buffer.len() == len);
				for &v in buffer.iter() { assert_eq!(v, buffer[0]); }
				buffer[0]
			})).try_into().unwrap():Box<_>),
			//eval(&jacobian, |buffer:&Buffer| buffer.map(device).unwrap()[0])
		);*/
		let all_same = |buffer:&Buffer| {
			let buffer = buffer.map(device).unwrap();
			assert!(buffer.len() == len);
			for &v in buffer.iter() { assert_eq!(v, buffer[0]); }
			buffer[0]
		};
		let gpu_transport = Transport{
			viscosity: all_same(&viscosity),
			thermal_conductivity: all_same(&thermal_conductivity),
			mixture_averaged_thermal_diffusion_coefficients: eval(&mixture_averaged_thermal_diffusion_coefficients, |buffer| all_same(buffer)),
		};
		if gpu_transport != transport { println!("{:?}\n{:?}", transport, gpu_transport); }
		println!("{:.0}K in {:.1}ms = {:.0}M/s", len as f32/1e3, time*1e3, (len as f32)/1e6/time);
	}
}
