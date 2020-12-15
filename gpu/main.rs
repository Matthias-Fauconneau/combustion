#![feature(default_free_fn,in_band_lifetimes, array_map, type_ascription,associated_type_bounds)]
mod vulkan;

#[fehler::throws(anyhow::Error)] fn main() {
	use {std::convert::TryInto, iter::{array_from_iter as from_iter, vec::{eval/*, generate*/}, box_collect}};
	let system = std::fs::read("CH4+O2.ron")?;
	const S : usize = 21;
	type Simulation<'t> = combustion::Simulation::<'t, S>;
	let Simulation{system, pressure_r, state: combustion::State{temperature, amounts}, ..} = Simulation::new(&system)?;
	let state = {use iter::into::IntoChain; from_iter([temperature].chain(amounts))};
	let f = system.dt_J(pressure_r, &state);
	let f = (f.0, /*f.1.concat().try_into().unwrap()*/);

	use vulkan::{Device, Buffer};
	let ref device = Device::new()?;
	let stride = 32;
	let len = 100000/stride*stride;
	let temperature = Buffer::new(device, (0..len).map(|_| temperature))?;
	let amounts_buffers = eval(amounts, |n| Buffer::new(device, (0..len).map(|_| n)).unwrap());
	let ref_amounts = box_collect(amounts_buffers.iter().map(|n| n));
	let d_temperature = Buffer::new(device, (0..len).map(|_| f64::NAN))?;
	let d_amounts = eval(amounts, |_| Buffer::new(device, (0..len).map(|_| f64::NAN)).unwrap());
	let ref_d_amounts = box_collect(d_amounts.iter().map(|n| n));
	let ref buffers = [&[&temperature] as &[_], &ref_amounts, &[&d_temperature] as &[_], &ref_d_amounts/*, &box_collect(jacobian.iter_mut().map(|n| n))*/];
	//let mut jacobian: [_;(2+S-1)*(2+S-1)] = generate(|_| Buffer::new(device, (0..len).map(|_| f64::NAN)).unwrap());
	let ref constants = [pressure_r];
	let ref pipeline = device.pipeline(constants, buffers)?;
	device.bind(pipeline.descriptor_set, buffers)?;
	let command_buffer  = device.command_buffer(pipeline, constants, stride, len)?;
	//let _warmup = device.submit_and_wait(command_buffer)?;
	for _ in 0..10 {
		let start = std::time::Instant::now();
		let time = device.submit_and_wait(command_buffer)?;
		let end = std::time::Instant::now();
		let cpu_time = (end-start).as_secs_f32();
		assert!(cpu_time > time);
		let gpu_f = (
			*(box_collect(std::iter::once(&d_temperature).chain(d_amounts.iter()).map(|buffer:&Buffer| {
				let buffer = buffer.map(device).unwrap();
				assert!(buffer.len() == len);
				for &v in buffer.iter() { assert_eq!(v, buffer[0]); }
				buffer[0]
			})).try_into().unwrap():Box<_>),
			//eval(&jacobian, |buffer:&Buffer| buffer.map(device).unwrap()[0])
		);
		use itertools::Itertools;
		if gpu_f != f { println!("{:.8e}\n{:.8e}", f.0.iter().format(" "), gpu_f.0.iter().format(" ")); }
		println!("{:.0}K in {:.1}ms = {:.0}M/s", len as f32/1e3, time*1e3, (len as f32)/1e6/time);
	}
}
