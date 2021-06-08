#![feature(bool_to_option)]#![allow(non_snake_case)]
use {fehler::throws, anyhow::Error};

use reaction::Simulation;
#[throws] fn main() {
	let model = std::fs::read("LiDryer.ron")?;
	let width = 32;
	let simulation = Simulation::new(&model, 512*width)?;
	let states_len = simulation.states_len();
	let Simulation{species_names, function, states, rates, pressure_Pa_R, temperature, mass_rate, energy_rate_R} = simulation;
	let values_count = function.dfg.values().count();
	let store_values = false;
	let function = cuda::cu(function, store_values)?;
	if std::fs::metadata("/var/tmp/main.cu").map_or(true, |cu|
		std::fs::read("/var/tmp/main.cu").unwrap() != function.as_bytes() ||
		std::fs::metadata("/var/tmp/main.ptx").map_or(true, |ptx| ptx.modified().unwrap() < cu.modified().unwrap())
	) {
		std::fs::write("/var/tmp/main.cu", &function)?;
		std::process::Command::new("nvcc").args(&["--ptx","/var/tmp/main.cu","-o","/var/tmp/main.ptx"]).spawn()?.wait()?.success().then_some(()).unwrap();
	}
	rustacuda::init(CudaFlags::empty()).unwrap();
	use rustacuda::{prelude::*, launch};
	let device = Device::get_device(0).unwrap();
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).unwrap();
	let module = Module::load_from_string(&std::ffi::CString::new(std::fs::read("/var/tmp/main.ptx").unwrap()).unwrap()).unwrap();
	let stream = Stream::new(/*irrelevant*/StreamFlags::NON_BLOCKING, None).unwrap();
	let mut device_states = DeviceBuffer::from_slice(&states).unwrap();
	let mut device_rates = DeviceBuffer::from_slice(&rates).unwrap();
	let start = std::time::Instant::now();
	if store_values {
		let values = vec![f64::NAN; values_count];
		let mut device_values = DeviceBuffer::from_slice(&values).unwrap();
		unsafe{launch!(module._Z6kerneldmmmmdddPd<<</*workgroupCount*/(states_len/width) as u32,/*workgroupSize*/width as u32, 0, stream>>>(
			pressure_Pa_R,
			device_states.as_device_ptr(),
			device_states.as_device_ptr().add(states_len),
			device_rates.as_device_ptr().add(states_len),
			device_rates.as_device_ptr(),
			temperature,
			1./mass_rate,
			1./energy_rate_R,
			device_values.as_device_ptr()
		)).unwrap()}
		stream.synchronize().unwrap();
		let mut values = values;
		device_values.copy_to(&mut values).unwrap();
		pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		std::fs::write("/var/tmp/values", as_bytes(&values))?;
	} else {
		unsafe{launch!(module._Z6kerneldmmmmddd<<</*workgroupCount*/(states_len/width) as u32,/*workgroupSize*/width as u32, 0, stream>>>(
			pressure_Pa_R,
			device_states.as_device_ptr(),
			device_states.as_device_ptr().add(states_len),
			device_rates.as_device_ptr().add(states_len),
			device_rates.as_device_ptr(),
			temperature,
			1./mass_rate,
			1./energy_rate_R
		)).unwrap()}
		stream.synchronize().unwrap();
	}
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f64();
	println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e3, (states_len as f64)/1e3/time);
	let mut rates = rates;
	device_rates.copy_to(&mut rates).unwrap();
	reaction::report(&species_names, &rates);
	println!("{:8} {:e}", "HRR", rates[0]);
}
