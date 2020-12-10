#![feature(default_free_fn,in_band_lifetimes, array_map, type_ascription)]
mod vulkan;

#[fehler::throws(anyhow::Error)] fn main() {
	use {std::convert::TryInto, iter::{array_from_iter as from_iter, vec::{eval, generate}, box_collect}};
	let system = std::fs::read("H2+O2.ron")?;
	const S : usize = 9;
	type Simulation<'t> = combustion::Simulation::<'t, S>;
	let Simulation{system, pressure_r, volume, state: combustion::State{temperature, amounts}, ..} = Simulation::new(&system)?;
	let state = {use iter::into::IntoChain; from_iter([temperature,volume].chain(amounts))};
	let f = system.dt(pressure_r, &state);
	let f = (f.0, f.1.concat().try_into().unwrap());

	use vulkan::{Device, Buffer};
	let ref device = Device::new()?;
	let stride = 1;
	let len = 1*stride;
	let [pressure_r, mut temperature, mut volume] = [pressure_r, temperature, volume].map(|initial| Buffer::new(device, (0..len).map(|_| initial)).unwrap());
	let mut amounts = eval(amounts, |n| Buffer::new(device, (0..len).map(|_| n)).unwrap());
	let mut jacobian: [_;(2+S-1)*(2+S-1)] = generate(|_| Buffer::new(device, (0..len).map(|_| f64::NAN)).unwrap());
	device.submit_and_wait(&[&[&pressure_r],&[&mut temperature],&[&mut volume], &box_collect(amounts.iter_mut().map(|n| n)), &box_collect(jacobian.iter_mut().map(|n| n))], stride, len)?;
	let gpu_f = (
		*(box_collect([temperature, volume].iter().chain(amounts.iter()).map(|buffer:&Buffer| buffer.map(device).unwrap()[0])).try_into().unwrap():Box<_>),
		eval(&jacobian, |buffer:&Buffer| buffer.map(device).unwrap()[0])
	);
	use itertools::Itertools;
	if gpu_f != f { println!("{:.8e}\n{:.8e}", f.0.iter().format(" "), gpu_f.0.iter().format(" ")); }
}
