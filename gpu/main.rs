#![feature(default_free_fn,in_band_lifetimes, array_map, type_ascription,associated_type_bounds)]
mod vulkan;

#[fehler::throws(anyhow::Error)] fn main() {
	use {std::convert::TryInto, iter::{array_from_iter as from_iter, vec::{eval/*, generate*/}, box_collect}};
	let system = std::fs::read("H2+O2.ron")?;
	const S : usize = 9;
	type Simulation<'t> = combustion::Simulation::<'t, S>;
	let Simulation{system, pressure_r, state: combustion::State{temperature, amounts}, ..} = Simulation::new(&system)?;
	let state = {use iter::into::IntoChain; from_iter([temperature].chain(amounts))};
	let f = system.dt_J(pressure_r, &state);
	let f = (f.0, /*f.1.concat().try_into().unwrap()*/);

	use vulkan::{Device, Buffer, BufferUsageFlags};
	let ref device = Device::new()?;
	let stride = 1792;
	let len = 65535*stride;
	let mut temperature = Buffer::new(device, BufferUsageFlags::STORAGE_BUFFER, (0..len).map(|_| temperature))?;
	let mut amounts = eval(amounts, |n| Buffer::new(device, BufferUsageFlags::STORAGE_BUFFER, (0..len).map(|_| n)).unwrap());
	//let mut jacobian: [_;(2+S-1)*(2+S-1)] = generate(|_| Buffer::new(device, (0..len).map(|_| f64::NAN)).unwrap());
	let start = std::time::Instant::now();
	let _time = device.submit_and_wait(&[pressure_r], &[&[&mut temperature], &box_collect(amounts.iter_mut().map(|n| n))/*, &box_collect(jacobian.iter_mut().map(|n| n))*/], stride, len)?;
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f32();
	let gpu_f = (
		*(box_collect([temperature].iter().chain(amounts.iter()).map(|buffer:&Buffer| buffer.map(device).unwrap()[0])).try_into().unwrap():Box<_>),
		//eval(&jacobian, |buffer:&Buffer| buffer.map(device).unwrap()[0])
	);
	use itertools::Itertools;
	if gpu_f != f { println!("{:.8e}\n{:.8e}", f.0.iter().format(" "), gpu_f.0.iter().format(" ")); }
	println!("{:.0}M in {:.0}ms = {:.0}M/s", len as f32/1e6, time*1e3, (len as f32)/1e6/time);
}
