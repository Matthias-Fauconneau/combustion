#![feature(type_ascription, array_map, array_methods, bindings_after_at, try_blocks)] #![allow(non_snake_case)]
#[fehler::throws(anyhow::Error)] fn main() {
	use iter::{box_collect, vec::ConstRange, into::map, array_from_iter as from_iter};
	let system = std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let Simulation{system, pressure_R, state: state@combustion::State{temperature, amounts}, ..} = Simulation::<35>::new(&system)?;
	let transport =  system.transport(pressure_R, &state);

	let _ : anyhow::Result<_> = try { std::process::Command::new("gpu-on").spawn()?.wait()? };
	//assert!(std::fs::read("/sys/devices/virtual/hwmon/hwmon4/temp9_input").is_ok());
	std::process::Command::new("nvidia-modprobe").spawn()?.wait()?;
	assert!(std::str::from_utf8(&std::fs::read("/proc/modules")?)?.lines().any(|line| line.starts_with("nvidia ")));
	use rustacuda::{prelude::*, launch};
	rustacuda::init(CudaFlags::empty())?;
	let device = Device::get_device(0)?;
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device)?;
	let module = Module::load_from_string(&std::ffi::CString::new(include_str!(concat!(env!("OUT_DIR"), "/main.ptx")))?)?;
	let stream = Stream::new(StreamFlags::NON_BLOCKING, None)?;

	let stride = 1;//256;
	let len = (1/*00_000*//stride)*stride;
	let mut temperature = DeviceBuffer::from_slice(&vec![temperature; len])?;
	let mut amounts_buffer = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|&n| std::iter::repeat(n).take(len)).flatten()))?;
	let mut d_temperature = DeviceBuffer::from_slice(&vec![f64::NAN; len])?;
	let mut d_amounts = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|_| std::iter::repeat(f64::NAN).take(len)).flatten()))?;
	let mut viscosity = DeviceBuffer::from_slice(&vec![f64::NAN; len])?;
	let mut thermal_conductivity = DeviceBuffer::from_slice(&vec![f64::NAN; len])?;
	let mut mixture_averaged_thermal_diffusion_coefficients = DeviceBuffer::from_slice(&box_collect(amounts.iter().map(|_| std::iter::repeat(f64::NAN).take(len)).flatten()))?;

	for _ in 0..1/*0*/ {
		let start = std::time::Instant::now();
		unsafe {
			launch!(module.rates_transport<<</*workgroupCount*/(len/stride) as u32,/*workgroupSize*/stride as u32, 0, stream>>>(len, pressure_R,
									 temperature.as_device_ptr(), amounts_buffer.as_device_ptr(), d_temperature.as_device_ptr(), d_amounts.as_device_ptr(), viscosity.as_device_ptr(), thermal_conductivity.as_device_ptr(), mixture_averaged_thermal_diffusion_coefficients.as_device_ptr()))?;
		}
		stream.synchronize()?;
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
				from_iter(map(ConstRange::<35>, move |n| host[n*len]))
			}
		};
		use combustion::Error;
		if transport.error(&gpu_transport) > 3e-6 { println!("{:?}\n{:?}", transport, gpu_transport); }
		println!("{:.0}K in {:.1}ms = {:.0}M/s", len as f32/1e3, time*1e3, (len as f32)/time/1e6);
	}
}
