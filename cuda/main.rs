#![feature(type_ascription, array_map, array_methods, bindings_after_at, try_blocks, unboxed_closures)] #![allow(non_snake_case)]

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
	let ref state = Simulation::new(&model)?.state;
	//#[cfg(feature="transport")] {
	let (_species_names, species) = Species::new(model.species);
	let ref transport_polynomials = species.transport_polynomials();
	let transport = benchmark!(transport::transport(&species.molar_mass, transport_polynomials, state), 1);

	let _ : std::io::Result<_> = try { std::process::Command::new("gpu-on").spawn()?.wait()? };
	let _ : std::io::Result<_> = try { std::process::Command::new("nvidia-modprobe").spawn()?.wait()? };
	if let Ok(modules) = std::fs::read("/proc/modules") { assert!(std::str::from_utf8(&modules).unwrap().lines().any(|line| line.starts_with("nvidia "))) };
	use rustacuda::{prelude::*, launch};
	rustacuda::init(CudaFlags::empty()).expect("CUDA");
	let device = Device::get_device(0).expect("device");
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).expect("context");
	let module = Module::load_from_string(&std::ffi::CString::new(include_str!(concat!(env!("OUT_DIR"), "/main.ptx"))).unwrap()).expect("module");
	let stream = Stream::new(StreamFlags::NON_BLOCKING, None).expect("stream");

	let stride = 1;//256;
	let len = (1/*00_000*//stride)*stride;
	if len < 1000 { println!("{}", len) } else { println!("{}K", len/1000); }
	let State{temperature, pressure, amounts, ..} = state;
	let mut temperature = DeviceBuffer::from_slice(&vec![*temperature; len]).unwrap();
	let mut amounts_buffer = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|&n| std::iter::repeat(n).take(len)).flatten())).unwrap();
	let mut d_temperature = DeviceBuffer::from_slice(&vec![f64::NAN; len]).unwrap();
	let mut d_amounts = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|_| std::iter::repeat(f64::NAN).take(len)).flatten())).unwrap();
	let mut viscosity = DeviceBuffer::from_slice(&vec![f64::NAN; len]).unwrap();
	let mut thermal_conductivity = DeviceBuffer::from_slice(&vec![f64::NAN; len]).unwrap();
	let mut mixture_averaged_thermal_diffusion_coefficients = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|_| std::iter::repeat(f64::NAN).take(len)).flatten())).unwrap();

	for _ in 0..1/*0*/ {
		let start = std::time::Instant::now();
		unsafe {
			launch!(module.rates_transport<<</*workgroupCount*/(len/stride) as u32,/*workgroupSize*/stride as u32, 0, stream>>>(len, *pressure,
									 temperature.as_device_ptr(), amounts_buffer.as_device_ptr(), d_temperature.as_device_ptr(), d_amounts.as_device_ptr(), viscosity.as_device_ptr(), thermal_conductivity.as_device_ptr(), mixture_averaged_thermal_diffusion_coefficients.as_device_ptr())).expect("launch");
		}
		stream.synchronize().expect("synchronize");
		let end = std::time::Instant::now();
		let time = (end-start).as_secs_f32();
		/*let gpu_f: [_;1+S-1] = from_iter([{
			let ref mut host = vec![0.; len];
			d_temperature.copy_to(host).unwrap();
			host[0]
		}].chain({
			let mut host = vec![0.; (S-1)*len];
			d_amounts.copy_to(&mut host).unwrap();
			generate(move |n| host[n*len])
		}));*/
		let gpu_transport = Transport{
			viscosity: {
				let ref mut host = vec![0.; len];
				viscosity.copy_to(host).unwrap();
				host[0]
			},
			thermal_conductivity: {
				let ref mut host = vec![0.; len];
				thermal_conductivity.copy_to(host).unwrap();
				host[0]
			},
			mixture_averaged_thermal_diffusion_coefficients: {
				let mut host = vec![0.; mixture_averaged_thermal_diffusion_coefficients.len()];
				mixture_averaged_thermal_diffusion_coefficients.copy_to(&mut host).unwrap();
				(0..len).map(|n| host[n*len]).collect()
			}
		};
		use AbsError;
		//if dbg!(transport.viscosity.error(&all_same(&gpu_transport.viscosity))) > 3e-6 { println!("{:?}\n{:?}", transport.viscosity, all_same(&viscosity)); }
		if transport.error(&gpu_transport) > 3e-6 { println!("{:?}\n{:?}", transport, gpu_transport); }
		println!("{:.1}ms\t{:.0}M/s", time*1e3, (len as f32)/time/1e6);
	}
}
