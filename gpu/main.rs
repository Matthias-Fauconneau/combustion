#![feature(default_free_fn, in_band_lifetimes, bindings_after_at, unboxed_closures)] #![allow(non_snake_case)]
mod vulkan;

fn time<T:Fn<()>>(task: T) -> (T::Output, f32) {
	let start = std::time::Instant::now();
	let result = task();
	(result, (std::time::Instant::now()-start).as_secs_f32())
}

fn print_time<T:Fn<()>>(task: T, header: impl std::fmt::Display) -> T::Output {
	let (result, time) = time(task);
	println!("{}: {:.1}ms", header, time*1e3);
	result
}

macro_rules! time { ($task:expr) => { print_time(|| { $task }, stringify!($task)) } }

#[fehler::throws(Box<dyn std::error::Error>)] fn main() {
	let model = &std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let model = model::Model::new(&model)?;
	let (species_names, species) = combustion::Species::new(&model.species);
	use reaction::*;
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(&species_names, r));
	//let width = 32;
	let states_len = 1; //((512*32)/width)*width;
	let instructions = rate::<_,{Property::Volume}>(&species, &*reactions, states_len)?;
	let main = std::str::from_utf8(&std::fs::read("main.comp")?)?.replace("#include \"instructions\"", &instructions).replace("float", "double");
	std::fs::write("/var/tmp/main.comp", &main)?;
  let main = shaderc::Compiler::new().unwrap().compile_into_spirv(&main, shaderc::ShaderKind::Compute, "main.comp", "main", None)?;

	use vulkan::{Device, Buffer};
	let ref device = Device::new()?;
	let stride = 1;
	let states_len = 1/stride*stride;
	let ref state = initial_state(&model);
	//let constant = state.constant::<{Property::Volume}>();
	let state : StateVector<{Property::Volume}> = state.into();
	let state = Buffer::new(device, iter::box_collect(state.iter().map(|&n| std::iter::repeat(n as f64/*f32*/).take(states_len)).flatten()).iter().copied())?;
	let HRR_dtω = Buffer::new(device, iter::box_collect(std::iter::repeat(f64/*f32*/::NAN).take((1+species.len()-1)*states_len)).iter().copied())?;
	let ref constants = [];
	let ref buffers = [&[&state] as &[_], &[&HRR_dtω] as &[_]];

	let ref pipeline = time!(device.pipeline(main.as_binary(), constants, buffers))?; // Compiles SPIRV -> Gen
	device.bind(pipeline.descriptor_set, buffers)?;
	let command_buffer  = device.command_buffer(pipeline, constants, stride, states_len)?;
	for _ in 0..1 {
		dbg!();
		let time = device.submit_and_wait(command_buffer)?;
		dbg!();
		let _all_same = |buffer:&Buffer| -> f64 {
			let buffer = buffer.map(device).unwrap();
			assert!(buffer.len() == states_len);
			for &v in buffer.iter() { assert_eq!(v, buffer[0]); }
			buffer[0]
		};
		println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e3, (states_len as f32)/1e3/time);
	}
}
