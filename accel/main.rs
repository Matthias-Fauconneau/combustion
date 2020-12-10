#![feature(type_ascription, array_map, array_methods)]
mod cuda;

pub const S : usize = 9;
pub const N : usize = 2+S-1;

mod kernel;

#[fehler::throws(anyhow::Error)] fn main() {
	use {std::convert::TryInto, iter::{array_from_iter as from_iter, vec::eval, box_collect}};
	let system = std::fs::read("H2+O2.ron")?;
	type Simulation<'t> = combustion::Simulation::<'t, S>;
	let Simulation{system, pressure_r, volume, state: combustion::State{temperature, amounts}, ..} = Simulation::new(&system)?;
	let state = {use iter::into::IntoChain; from_iter([temperature,volume].chain(amounts))};
	let f = system.dt(pressure_r, &state);

	use accel::{DeviceMemory, Allocatable};
	let ref device = cuda::Device::new()?;
	let stride = 1;
	let len = 1*stride;
	let [pressure_r, mut temperature, mut volume] = [pressure_r, temperature, volume].map(|initial| DeviceMemory::<f32>::from_elem(device, len, initial as f32));
	let mut amounts = eval(amounts, |n| DeviceMemory::<f32>::from_elem(device, len, n as f32));
	//let mut jacobian: [[_;N]; N] = generate(|_| generate(|_| DeviceMemory::<f32>::from_elem(device, len, f32::NAN)));
	kernel::dt(device, 1, len, (len,
		pressure_r.as_ptr(), temperature.as_mut_ptr(), volume.as_mut_ptr(), box_collect(amounts.iter_mut().map(|n| n.as_mut_ptr())).as_ref(),//.as_slice(),
		//box_collect(jacobian.iter_mut().map(|j| box_collect(j.map(|j| j.as_mut_ptr())).as_slice())).as_slice()
	))?;
	let gpu_f : [_;2+S-1] = *(box_collect([temperature, volume].iter().chain(amounts.iter()).map(|buffer| buffer[0] as f64)).try_into().unwrap():Box<_>);//, eval(&jacobian, |row| eval(row, |buffer| buffer[0] as f64)))
	use itertools::Itertools;
	if gpu_f != f.0 { println!("{:.8e}\n{:.8e}", f.0.iter().format(" "), gpu_f.iter().format(" ")); }
}
