#![feature(type_ascription, array_map, array_methods)]
#[fehler::throws(anyhow::Error)] fn main() {
	use {std::convert::TryInto, iter::{array_from_iter as from_iter, box_collect}};
	let system = std::fs::read("H2+O2.ron")?;
	pub const S : usize = 9;

	type Simulation<'t> = combustion::Simulation::<'t, S>;
	let Simulation{system, pressure_r, state: combustion::State{temperature, amounts}, ..} = Simulation::new(&system)?;
	let state = {use iter::into::IntoChain; from_iter([temperature].chain(amounts))};
	let f = system.dt_J(pressure_r, &state);

	std::process::Command::new("gpu-on").spawn()?.wait()?;
	assert!(std::fs::read("/sys/devices/virtual/hwmon/hwmon4/temp9_input").is_ok());
	std::process::Command::new("nvidia-modprobe").spawn()?.wait()?;
	assert!(std::str::from_utf8(&std::fs::read("/proc/modules")?)?.lines().any(|line| line.starts_with("nvidia ")));
	use rustacuda::{prelude::*, launch};
	rustacuda::init(CudaFlags::empty())?;
	let device = Device::get_device(0)?;
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device)?;
	let module = Module::load_from_string(&std::ffi::CString::new(include_str!(concat!(env!("OUT_DIR"), "/main.ptx")))?)?;
	let stream = Stream::new(StreamFlags::NON_BLOCKING, None)?;

	let stride = 1;
	let len = 1*stride;
	let mut temperature = DeviceBuffer::from_slice(&vec![temperature; len])?;
	let mut amounts_buffer = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|&n| std::iter::repeat(n).take(len)).flatten()))?;
	let mut d_temperature = DeviceBuffer::from_slice(&vec![f64::NAN; len])?;
	let mut d_amounts = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|_| std::iter::repeat(f64::NAN).take(len)).flatten()))?;
	unsafe { launch!(module.dt<<<1, 1, 0, stream>>>(len, pressure_r, temperature.as_device_ptr(), amounts_buffer.as_device_ptr(), d_temperature.as_device_ptr(), d_amounts.as_device_ptr()))?; }
	stream.synchronize()?;
	let gpu_f : [_;1+S-1] = *(box_collect(std::iter::once({let ref mut host = [0.; 1]; d_temperature.copy_to(host).unwrap(); host[0] }).chain({
		let mut host = vec![0.; (S-1)*len]; d_amounts.copy_to(&mut host).unwrap(); (0..(S-1)).map(move |n| host[n*len])
	})).try_into().unwrap():Box<_>);
	use itertools::Itertools;
	if gpu_f != f.0 { println!("{:.8e}\n{:.8e}", f.0.iter().format(" "), gpu_f.iter().format(" ")); }
}
