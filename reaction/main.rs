#![allow(non_snake_case)]#![feature(bool_to_option,assert_matches)]

pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
pub fn as_bytes_mut<T>(slice: &mut [T]) -> &mut [u8] { unsafe{std::slice::from_raw_parts_mut(slice.as_mut_ptr() as *mut u8, slice.len() * std::mem::size_of::<T>())} }

use cranelift_codegen::{ir::{function::Function, immediates::Ieee64}, data_value::DataValue};

enum Argument<'t> { Value(DataValue), Ref(&'t [u8]), Mut(&'t mut[u8]) }

fn interpret(function: &Function, arguments: &mut [Argument]) {
	use cranelift_interpreter::interpreter::{InterpreterState, Interpreter};
	use Argument::*;
	let mut state = InterpreterState::default();
	let arguments_values = iter::map(arguments.iter(), |v| {
		match v {
			Value(v) => v.clone(),
			Ref(v) => { let base = state.heap.len(); state.heap.extend(v.iter()); DataValue::I32(base as i32) }
			Mut(v) => { let base = state.heap.len(); state.heap.extend(v.iter()); DataValue::I32(base as i32) }
		}
	});
	let mut interpreter = Interpreter::new(state);
	assert!(interpreter.call(function, &arguments_values).unwrap().unwrap_return() == &[]);
	let mut base = 0;
	let state = interpreter.state;
	for v in arguments.iter_mut() {
		match v {
			Value(_) => (),
			Ref(v) => base += v.len(),
			Mut(v) => { v.copy_from_slice(&state.heap[base..base+v.len()]); base += v.len(); }
		}
	}
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = &std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let model = model::Model::new(&model)?;
	let (species_names, species) = combustion::Species::new(&model.species);
	use reaction::*;
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(&species_names, r));
	let width = 32;
	let states_len = ((1*32)/width)*width;
	let function = rate::<_,{Property::Pressure}>(&species, &*reactions, states_len)?;

	//std::fs::write("/var/tmp/main.cu", &function)?;
	//std::process::Command::new("nvcc").args(&["--ptx","/var/tmp/main.cu","-o","/var/tmp/main.ptx"]).spawn()?.wait()?.success().then_some(()).unwrap();

	let ref reference_state = initial_state(&model);
	let length = 1f64;
	let velocity = 1f64;
	let time = length / velocity;
	let total_concentration = reference_state.pressure_R / reference_state.temperature;
	let total_amount = reference_state.amounts.iter().sum::<f64>();
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
	let total_amount = state.amounts.iter().sum::<f64>();
	let mole_fractions = iter::map(&*state.amounts, |n| n / total_amount);
	let molar_mass = dot(&*mole_fractions, &*species.molar_mass);
	let mass_fractions = iter::map(mole_fractions.iter().zip(species.molar_mass.iter()), |(x,m)| x * m / molar_mass);
	let ref state = [&[state.temperature / reference_state.temperature] as &[_], &*mass_fractions].concat();

	/*rustacuda::init(CudaFlags::empty()).unwrap();
	use rustacuda::{prelude::*, launch};
	let device = Device::get_device(0).unwrap();
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).unwrap();
	let module = Module::load_from_string(&std::ffi::CString::new(std::fs::read("/var/tmp/main.ptx").unwrap()).unwrap()).unwrap();
	let stream = Stream::new(StreamFlags::NON_BLOCKING, None).unwrap();*/

	let states = iter::box_collect(state.iter().map(|&s| std::iter::repeat(s as f64).take(states_len)).flatten());
	let mut rates = vec![f64::NAN; (1/*2*/+species.len())*states_len];

	for _ in 0..1 {
		let start = std::time::Instant::now();
		if false {
			/*let mut device_states = DeviceBuffer::from_slice(&states).unwrap();
			let mut device_rates = DeviceBuffer::from_slice(&rates).unwrap();
			unsafe{launch!(module._Z6kerneldPdS_S_ddS_<<</*workgroupCount*/(states_len/width) as u32,/*workgroupSize*/width as u32, 0, stream>>>(
				pressure_R,
				device_states.as_device_ptr(),
				device_states.as_device_ptr().add(states_len),
				device_rates.as_device_ptr().add(states_len),
				mass_production_rate_factor,
				heat_release_rate_factor,
				device_rates.as_device_ptr()
			)).unwrap()}
			device_rates.copy_to(&mut rates).unwrap();
			stream.synchronize().unwrap();*/
		} else {
			let (temperature, mass_fractions) = states.split_at(states_len);
			let (heat_release_rate, mass_production_rates) = rates.split_at_mut(states_len);
			use Argument::*;
			use DataValue::*;
			interpret(&function, &mut [
				Value(I32(0)),
				Value(F64(Ieee64::with_float(pressure_R))),
				Ref(as_bytes(temperature)),
				Ref(as_bytes(mass_fractions)),
				Mut(as_bytes_mut(mass_production_rates)),
				Value(F64(Ieee64::with_float(mass_production_rate_factor))),
				Value(F64(Ieee64::with_float(heat_release_rate_factor))),
				Mut(as_bytes_mut(heat_release_rate)),
			]);
		}
		let end = std::time::Instant::now();
		#[track_caller] fn all_same(slice: &[f64], stride: usize) -> Box<[f64]> {
      assert_eq!(slice.len()%stride, 0);
      iter::eval(slice.len()/stride, |i| {
        let slice = &slice[i*stride..(i+1)*stride];
        for &v in slice.iter() { assert_eq!(v, slice[0]); }
        slice[0]
      })
    }
		for ((mass_production_rate, name), molar_mass) in all_same(&rates, states_len).iter().skip(1).zip(species_names.iter()).zip(species.molar_mass.iter()) {
			let mass_production_rate = mass_production_rate * time/density * 1e8 / molar_mass;
			println!("{:5} {:e}", name, mass_production_rate)
		}
		let time = (end-start).as_secs_f64();
		println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e3, (states_len as f64)/1e3/time);
	}

	if false {
		use itertools::Itertools;
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
