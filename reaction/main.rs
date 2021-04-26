//#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked, bool_to_option, non_ascii_idents, array_methods, slice_pattern)]
//#![feature(type_ascription, array_map, array_methods, bindings_after_at, try_blocks, unboxed_closures)] #![allow(non_snake_case)]
#![feature(unboxed_closures, bool_to_option)]
//use std::{ops::Deref, convert::TryInto};

/*fn time<T:Fn<()>>(task: T) -> (T::Output, f32) {
	let start = std::time::Instant::now();
	let result = task();
	(result, (std::time::Instant::now()-start).as_secs_f32())
}

fn benchmark<T>(task: impl Fn()->T, times: usize, header: impl std::fmt::Display) -> T {
	let result = task();
	println!("{}: {:.1}ms", header, time(|| for _ in 0..times { task(); }).1/(times as f32)*1e3);
	result
}*/

//macro_rules! benchmark { ($task:expr, $times:expr) => { benchmark(|| { $task }, $times, stringify!($task)) } }

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = &std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let model = model::Model::new(&model)?;
	let (species_names, species) = combustion::Species::new(&model.species);
	use reaction::*;
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(&species_names, r));
	let width = 1;
	let states_len = ((1)/width)*width;
	let instructions = rate::<_,{Property::Pressure}>(&species, &*reactions, states_len)?;

	use itertools::Itertools;
	std::fs::write("/var/tmp/main.cu", &instructions)?;
	std::process::Command::new("nvcc").args(&["--ptx","/var/tmp/main.cu","-o","/var/tmp/main.ptx"]).spawn()?.wait()?.success().then_some(()).unwrap();

	let ref reference_state = initial_state(&model);
	let length = 1.;
	let velocity = 1.;
	let time = length / velocity;
	let total_concentration = reference_state.pressure_R / reference_state.temperature;
	let total_amount = reference_state.amounts.iter().sum();
	let mole_fractions = iter::map(&*reference_state.amounts, |n| n / total_amount);
	let dot = |a:&[_],b:&[_]| iter::dot(a.iter().copied().zip(b.iter().copied()));
	let molar_mass = dot(&*mole_fractions, &*species.molar_mass);
	let density = total_concentration * molar_mass;
	let molar_heat_capacity_R = iter::dot(mole_fractions.iter().copied().zip(
		species.thermodynamics.iter().map(|specie| specie.molar_heat_capacity_at_constant_pressure_R(reference_state.temperature))
	));
	let mass_production_rate_factor = time / density;
	let heat_release_rate_factor = 1. / (molar_heat_capacity_R * (K*NA) * reference_state.temperature /* time / density*/);

	let ref state = reference_state;
	let pressure_R = state.pressure_R / reference_state.pressure_R;
	let total_amount = state.amounts.iter().sum();
	let mole_fractions = iter::map(&*state.amounts, |n| n / total_amount);
	let molar_mass = dot(&*mole_fractions, &*species.molar_mass);
	let mass_fractions = iter::map(mole_fractions.iter().zip(species.molar_mass.iter()), |(x,m)| x * m / molar_mass);
	let ref state = [&[state.temperature / reference_state.temperature] as &[_], &*mass_fractions].concat();

	rustacuda::init(CudaFlags::empty()).unwrap();
	use rustacuda::{prelude::*, launch, memory};
	let device = Device::get_device(0).unwrap();
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).unwrap();
	let module = Module::load_from_string(&std::ffi::CString::new(std::fs::read("/var/tmp/main.ptx").unwrap()).unwrap()).unwrap();
	let stream = Stream::new(StreamFlags::NON_BLOCKING, None).unwrap();

	let mut states = DeviceBuffer::from_slice(&iter::box_collect(state.iter().map(|&s| std::iter::repeat(s as f32).take(states_len)).flatten())).unwrap();
	let mut rates = DeviceBuffer::from_slice(&vec![f32::NAN; (1/*2*/+species.len())*states_len]).unwrap();

	module.get_global::<f64>(std::ffi::CStr::from_bytes_with_nul(b"p1\0")?)?.copy_from(&pressure_R)?;
	module.get_global::<*const f64>(std::ffi::CStr::from_bytes_with_nul(b"p2\0")?)?.copy_from(states.as_device_ptr())?;
	module.get_global::<*const f64>(std::ffi::CStr::from_bytes_with_nul(b"p3\0")?)?.copy_from(states.as_device_ptr()+states_len)?;
	module.get_global::<*mut f64>(std::ffi::CStr::from_bytes_with_nul(b"p4\0")?)?.copy_from(rates.as_device_ptr()+states_len)?;
	module.get_global::<f64>(std::ffi::CStr::from_bytes_with_nul(b"p5\0")?)?.copy_from(mass_production_rate_factor)?;
	module.get_global::<f64>(std::ffi::CStr::from_bytes_with_nul(b"p6\0")?)?.copy_from(heat_release_rate_factor)?;
	module.get_global::<*mut f64>(std::ffi::CStr::from_bytes_with_nul(b"p7\0")?)?.copy_from(rates.as_device_ptr())?;

	for _ in 0..1 {
		let start = std::time::Instant::now();
		unsafe{launch!(module.reaction<<</*workgroupCount*/(states_len/width) as u32,/*workgroupSize*/width as u32, 0, stream>>>()).unwrap()}
		stream.synchronize().unwrap();
		let end = std::time::Instant::now();
		let time = (end-start).as_secs_f32();
		#[track_caller] fn all_same(slice: &memory::DeviceSlice<f32>) -> f32 {
			let ref mut host = vec![0.; slice.len()];
			slice.copy_to(host).unwrap();
			for &v in host.iter() { assert_eq!(v, host[0]); }
			host[0]
		}

		/*pub trait Error { fn error(&self, o: &Self) -> f64; }
		impl Error for f64 { fn error(&self, o: &Self) -> f64 { num::relative_error(*self, *o) } }
		impl Error for [f64] { fn error(&self, o: &Self) -> f64 { self.iter().zip(o.iter()).map(|(s,o)| s.error(o)).reduce(f64::max).unwrap() } }
		fn check<T: std::fmt::Debug+Error+?Sized>(label: &str, reference: &T, value: &T, tolerance: f64) {
			let error = Error::error(reference, value);
			if error > tolerance { println!("{}\n{:?}\n{:?}\n{:e}", label, reference, value, error); }
		}
		check("Reaction", &cpu_dtT_dtS_dtn, &all_same(dtT_dtS_dtn), 0.);*/
		println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e3, (states_len as f32)/1e3/time);
	}

	if true {
		std::fs::write("/var/tmp/gri.h", [
			format!("const int nSpecies = {};\n", species.len()),
			format!("const char* species[] = {{\"{}\"}};\n", species_names.iter().format("\", \"",)),
			format!("const double molar_mass[] = {{ {} }};\n", species.molar_mass.iter().format(", ")),
			format!("const double temperature_split[] = {{ {} }};\n", species.thermodynamics.iter().map(|_| NASA7::temperature_split).format(", ")),
			format!("const double NASA7[][2][7] = {{ {} }};\n", species.thermodynamics.iter().map(|a| format!("{{ {} }}", a.0.iter().map(|a| format!("{{ {} }}", a.iter().format(", "))).format(", "))).format(",\n")),
			].concat())?;
	}

	Ok(())
}
