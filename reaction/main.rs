//#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked, bool_to_option, non_ascii_idents, array_methods, slice_pattern)]
//#![feature(type_ascription, array_map, array_methods, bindings_after_at, try_blocks, unboxed_closures)] #![allow(non_snake_case)]
#![feature(unboxed_closures, bool_to_option)]
//use std::{ops::Deref, convert::TryInto};

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

//macro_rules! benchmark { ($task:expr, $times:expr) => { benchmark(|| { $task }, $times, stringify!($task)) } }

use rustacuda::{prelude::*, launch, memory::DeviceSlice};

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = &std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let model = model::Model::new(&model)?;
	let (species_names, species) = combustion::Species::new(&model.species);
	use reaction::*;
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(&species_names, r));
	let width = 32;
	let states_len = ((512*32)/width)*width;
	let instructions = rate::<_,{Property::Volume}>(&species, &*reactions, states_len)?;

	std::fs::write("/var/tmp/instructions", instructions)?;
	std::process::Command::new("nvcc").args(&["-I/var/tmp", "--ptx","main.cu","-o","/var/tmp/main.ptx"]).spawn()?.wait()?.success().then_some(()).unwrap();

	let ref state = initial_state(&model);
	rustacuda::init(CudaFlags::empty()).unwrap();
	let device = Device::get_device(0).unwrap();
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).unwrap();
	let module = Module::load_from_string(&std::ffi::CString::new(std::fs::read("/var/tmp/main.ptx").unwrap()).unwrap()).unwrap();
	//module.get_global::<[f32; CONSTANTS_LEN]>(std::ffi::CStr::from_bytes_with_nul(b"constants\0")?)?.copy_from(constants.deref().try_into()?)?;
	let stream = Stream::new(StreamFlags::NON_BLOCKING, None).unwrap();

	let constant = state.constant();
	let state : StateVector<{Property::Volume}> = state.into();
	let mut state = DeviceBuffer::from_slice(&iter::box_collect(state.iter().map(|&n| std::iter::repeat(n as f32).take(states_len)).flatten())).unwrap();
	let mut dtT_dtS_dtn = DeviceBuffer::from_slice(&vec![f32::NAN; (2+species.len())*states_len]).unwrap();

	for _ in 0..1 {
		let start = std::time::Instant::now();
		unsafe {
			launch!(module.reaction<<</*workgroupCount*/(states_len/width) as u32,/*workgroupSize*/width as u32, 0, stream>>>(
				constant as f32, state.as_device_ptr(), dtT_dtS_dtn.as_device_ptr()
			)).unwrap()
		}
		stream.synchronize().unwrap();
		let end = std::time::Instant::now();
		let time = (end-start).as_secs_f32();
		#[track_caller] fn all_same(slice: &DeviceSlice<f32>) -> f32 {
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
		check("Reaction", &cpu_dtT_dtS_dtn, &dtT_dtS_dtn, 0.);*/
		println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e3, (states_len as f32)/1e3/time);
	}
	Ok(())
}
