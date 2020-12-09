#![feature(default_free_fn,in_band_lifetimes, array_map, type_ascription)]
mod vulkan; use vulkan::{Device, Buffer};

#[fehler::throws(anyhow::Error)] fn main() {
	use {std::convert::TryInto, iter::{array_from_iter as from_iter, vec::{eval, generate}, box_collect}};
	let system = std::fs::read("H2+O2.ron")?;
	type Simulation<'t> = combustion::Simulation::<'t, 9>;
	let Simulation{system, pressure, volume, state: combustion::State{temperature, amounts}, ..} = Simulation::new(&system)?;
	let state = {use iter::into::IntoChain; from_iter([temperature,volume].chain(amounts))};
	let f = system.dt(pressure, &state);
	let f = (f.0, f.1.concat().try_into().unwrap());

	let ref device = Device::new()?;
	let stride = 64;
	let len = 1*stride;
	let [mut pressure, mut temperature, mut volume] = [pressure, temperature, volume].map(|initial| Buffer::new(device, (0..len).map(|_| initial)).unwrap());
	let mut amounts = eval(amounts, |n| Buffer::new(device, (0..len).map(|_| n)).unwrap());
	let jacobian:[_;10*10] = generate(|_| Buffer::new(device, (0..len).map(|_| 0.)).unwrap());
	device.submit_and_wait(&[&[&mut pressure],&[&mut temperature],&[&mut volume], &box_collect(amounts.iter_mut().map(|n| n))], stride, len)?;
	let gpu_f = (
		*(box_collect([temperature, volume].iter().chain(amounts.iter()).map(|buffer:&Buffer| buffer.map(device).unwrap()[0])).try_into().unwrap():Box<_>),
		eval(&jacobian, |buffer:&Buffer| buffer.map(device).unwrap()[0])
	);
	use itertools::Itertools;
	if gpu_f != f { println!("{:.8e}\n{:.8e}", f.0.iter().format(" "), gpu_f.0.iter().format(" ")); }
}
