#![feature(type_ascription, array_map, array_methods, bindings_after_at, try_blocks, unboxed_closures)] #![allow(non_snake_case)]
use std::{ops::Deref, convert::TryInto};
mod parse;

fn time<T:Fn<()>>(task: T) -> (T::Output, f32) {
	let start = std::time::Instant::now();
	let result = task();
	(result, (std::time::Instant::now()-start).as_secs_f32())
}

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
	let ref state = initial_state(&model);
	//#[cfg(feature="transport")] {
	let (_species_names, species) = combustion::Species::new(&model.species);
	//dbg!(species.len());
	let transport_polynomials = species.transport_polynomials();
	assert_eq!(state.volume, 1.);
	let ref transport = benchmark!(transport::transport(&species.molar_mass, &transport_polynomials, state), 1);

	/*let _ : std::io::Result<_> = try { std::process::Command::new("gpu-on").spawn()?.wait()? };
	let _ : std::io::Result<_> = try { std::process::Command::new("nvidia-modprobe").spawn()?.wait()? };
	if let Ok(modules) = std::fs::read("/proc/modules") { assert!(std::str::from_utf8(&modules).unwrap().lines().any(|line| line.starts_with("nvidia "))) };*/
	use rustacuda::{prelude::*, launch, memory::DeviceSlice};
	rustacuda::init(CudaFlags::empty()).expect("CUDA");
	let device = Device::get_device(0).expect("device");
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).expect("context");
	let module = Module::load_from_string(&std::ffi::CString::new(include_str!(concat!(env!("OUT_DIR"), "/main.ptx"))).unwrap()).expect("module");
	const N: usize = parse::parse(env!("SPECIES"));
	use iter::map;
	let molar_mass = map(species.molar_mass.iter(), |&v| v as f32);
	let TransportPolynomials{sqrt_viscosity_T14, thermal_conductivity_T12, binary_thermal_diffusion_coefficients_T32} = &transport_polynomials;
	let sqrt_viscosity_T14 = map(sqrt_viscosity_T14.iter(), |e| e.map(|v| v as f32) );
	let thermal_conductivity_T12 = map(thermal_conductivity_T12.iter(), |e| e.map(|v| v as f32) );
	let binary_thermal_diffusion_coefficients_T32 = map(binary_thermal_diffusion_coefficients_T32.iter(), |e| map(e.iter(), |e| e.map(|v| v as f32) ) );
	module.get_global::<[f32; N]>(std::ffi::CStr::from_bytes_with_nul(b"molar_mass\0")?)?.copy_from(molar_mass.deref().try_into()?)?;
	module.get_global::<[[f32; 5]; N]>(std::ffi::CStr::from_bytes_with_nul(b"sqrt_viscosity_T14\0")?)?.copy_from(sqrt_viscosity_T14.deref().try_into()?)?;
	module.get_global::<[[f32; 5]; N]>(std::ffi::CStr::from_bytes_with_nul(b"thermal_conductivity_T12\0")?)?.copy_from(thermal_conductivity_T12.deref().try_into()?)?;
	module.get_global::<[[[f32; 5]; N]; N]>(std::ffi::CStr::from_bytes_with_nul(b"binary_thermal_diffusion_coefficients_T32\0")?)?
		.copy_from(&iter::from_iter_(binary_thermal_diffusion_coefficients_T32.iter().map(|e| (*e.deref()).try_into().unwrap())))?;

	let stream = Stream::new(StreamFlags::NON_BLOCKING, None).expect("stream");

	let stride = 32;
	let len = ((512*32)/stride)*stride;
	let State{temperature, pressure_R, amounts, ..} = state;
	let mut temperature = DeviceBuffer::from_slice(&vec![(*temperature) as f32; len]).unwrap();
	let mut amounts_buffer = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|&n| std::iter::repeat(n as f32).take(len)).flatten())).unwrap();
	let mut viscosity = DeviceBuffer::from_slice(&vec![f32::NAN; len]).unwrap();
	let mut thermal_conductivity = DeviceBuffer::from_slice(&vec![f32::NAN; len]).unwrap();
	let mut mixture_molar_averaged_thermal_diffusion_coefficients =
		DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|_| std::iter::repeat(f32::NAN).take(len)).flatten())).unwrap();

	for _ in 0..1 {
		let start = std::time::Instant::now();
		unsafe {
			launch!(module.rates_transport<<</*workgroupCount*/(len/stride) as u32,/*workgroupSize*/stride as u32, 0, stream>>>(
				len, (*pressure_R) as f32,
				temperature.as_device_ptr(), amounts_buffer.as_device_ptr(), viscosity.as_device_ptr(),
				thermal_conductivity.as_device_ptr(), mixture_molar_averaged_thermal_diffusion_coefficients.as_device_ptr())).expect("launch");
		}
		stream.synchronize().expect("synchronize");
		let end = std::time::Instant::now();
		let time = (end-start).as_secs_f32();
		#[track_caller] fn all_same(slice: &DeviceSlice<f32>) -> f32 {
			let ref mut host = vec![0.; slice.len()];
			slice.copy_to(host).unwrap();
			for &v in host.iter() { assert_eq!(v, host[0]); }
			host[0]
		}
		//map(slices, all_same)
		#[track_caller] fn map_all_same(slice: &DeviceSlice<f32>, stride: usize) -> Box<[f32]> {
			assert_eq!(slice.len()%stride, 0);
			let ref mut host = vec![0.; slice.len()];
			slice.copy_to(host).unwrap();
			iter::eval(slice.len()/stride, |i| {
				let host = &host[i*stride..(i+1)*stride];
				for &v in host.iter() { assert_eq!(v, host[0]); }
				host[0]
			})
		}

		let gpu_transport = Transport{
			viscosity: all_same(&viscosity) as f64,
			thermal_conductivity: all_same(&thermal_conductivity) as f64,
			mixture_molar_averaged_thermal_diffusion_coefficients: iter::map(map_all_same(&mixture_molar_averaged_thermal_diffusion_coefficients, len).iter(), |&v| v as f64)
		};

		pub trait Error { fn error(&self, o: &Self) -> f64; }
		impl Error for f64 { fn error(&self, o: &Self) -> f64 { num::relative_error(*self, *o) } }
		impl Error for [f64] { fn error(&self, o: &Self) -> f64 { self.iter().zip(o.iter()).map(|(s,o)| s.error(o)).reduce(f64::max).unwrap() } }
		fn check<T: std::fmt::Debug+Error+?Sized>(label: &str, reference: &T, value: &T, tolerance: f64) {
			let error = Error::error(reference, value);
			if error > tolerance { println!("{}\n{:?}\n{:?}\n{:e}", label, reference, value, error); }
		}
		check("Viscosity", &transport.viscosity, &gpu_transport.viscosity, 2e-6);
		check("Thermal Conductivity", &transport.thermal_conductivity, &gpu_transport.thermal_conductivity, 7e-6);
		check("Mixture molar averaged thermal diffusion coefficients", transport.mixture_molar_averaged_thermal_diffusion_coefficients.deref(),
																																															  gpu_transport.mixture_molar_averaged_thermal_diffusion_coefficients.deref(), 9e-6);
		println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", len as f32/1e3, time*1e3, time/(len as f32)*1e3, (len as f32)/1e3/time);
	}
}