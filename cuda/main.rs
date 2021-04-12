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
	//dbg!(species.len());
	let ref transport_polynomials = species.transport_polynomials();
	//transport_polynomials.binary_thermal_diffusion_coefficients_T32 // 5x53x(2+53)x4B/f32 = 56K (NV CUDA CC3 constant memory size: 64K)
	let transport = benchmark!(transport::transport(&species.molar_mass, transport_polynomials, state), 1);

	/*let _ : std::io::Result<_> = try { std::process::Command::new("gpu-on").spawn()?.wait()? };
	let _ : std::io::Result<_> = try { std::process::Command::new("nvidia-modprobe").spawn()?.wait()? };
	if let Ok(modules) = std::fs::read("/proc/modules") { assert!(std::str::from_utf8(&modules).unwrap().lines().any(|line| line.starts_with("nvidia "))) };*/
	use rustacuda::{prelude::*, launch, memory::DeviceSlice};
	rustacuda::init(CudaFlags::empty()).expect("CUDA");
	let device = Device::get_device(0).expect("device");
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).expect("context");
	println!("CUDA...");
	let module = Module::load_from_string(&std::ffi::CString::new(include_str!(concat!(env!("OUT_DIR"), "/main.ptx"))).unwrap()).expect("module");
	let stream = Stream::new(StreamFlags::NON_BLOCKING, None).expect("stream");

	let stride = 1;//256;
	let len = (1/*00_000*//stride)*stride;
	//if len < 1000 { println!("{}", len) } else { println!("{}K", len/1000); }
	let State{temperature, pressure, amounts, ..} = state;
	let mut temperature = DeviceBuffer::from_slice(&vec![*temperature as f32; len]).unwrap();
	let mut amounts_buffer = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|&n| std::iter::repeat(n as f32).take(len)).flatten())).unwrap();
	let mut d_temperature = DeviceBuffer::from_slice(&vec![f32::NAN; len]).unwrap();
	let mut d_amounts = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|_| std::iter::repeat(f32::NAN).take(len)).flatten())).unwrap();
	let mut viscosity = DeviceBuffer::from_slice(&vec![f32::NAN; len]).unwrap();
	let mut thermal_conductivity = DeviceBuffer::from_slice(&vec![f32::NAN; len]).unwrap();
	let mut mixture_averaged_thermal_diffusion_coefficients = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|_| std::iter::repeat(f32::NAN).take(len)).flatten())).unwrap();

	for _ in 0..1/*0*/ {
		let start = std::time::Instant::now();
		println!("CUDA...");
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
		#[track_caller] fn all_same(slice: &DeviceSlice<f32>) -> f32 {
			let ref mut host = vec![0.; slice.len()];
			slice.copy_to(host).unwrap();
			for &v in host.iter() { assert_eq!(v, host[0]); }
			host[0]
		}

		use RelError;

		if transport.viscosity.error(&(all_same(&viscosity) as f64)) > 2e-6 || true {
			println!("{:?}\n{:?}\n{:e}", transport.viscosity, all_same(&viscosity), transport.viscosity.error(&(all_same(&viscosity) as f64)));
		}

		if transport.thermal_conductivity.error(&(all_same(&thermal_conductivity) as f64)) > 3e-5 || true {
			println!("{:?}\n{:?}\n{:e}", transport.thermal_conductivity, all_same(&thermal_conductivity), transport.thermal_conductivity.error(&(all_same(&thermal_conductivity) as f64)));
		}

		let gpu_transport = Transport{
			viscosity: all_same(&viscosity) as f64,
			thermal_conductivity: all_same(&thermal_conductivity) as f64,
			mixture_averaged_thermal_diffusion_coefficients: //box_collect(mixture_averaged_thermal_diffusion_coefficients.iter().map(|buffer| all_same(buffer))),
				{
					let ref slice = mixture_averaged_thermal_diffusion_coefficients;
					let ref mut host = vec![0.; slice.len()];
					slice.copy_to(host).unwrap();
					(0..len).map(|n| host[n*len] as f64).collect()
				},
		};

		if dbg!(transport.error(&gpu_transport)) > 3e-5 {
			println!("{:?}\n{:?}\n{:e}", transport, gpu_transport, transport.error(&gpu_transport));
		}

		println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", len as f32/1e3, time*1e3, time/(len as f32)*1e3, (len as f32)/1e3/time);
	}
}
