#![allow(non_snake_case)]#![feature(bool_to_option,assert_matches,default_free_fn)]
use std::default::default;

pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
pub fn as_bytes_mut<T>(slice: &mut [T]) -> &mut [u8] { unsafe{std::slice::from_raw_parts_mut(slice.as_mut_ptr() as *mut u8, slice.len() * std::mem::size_of::<T>())} }
pub fn from_bytes<T>(slice: &[u8]) -> &[T] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const T, slice.len() / std::mem::size_of::<T>())} }

use cranelift_codegen::{ir::{function::Function, immediates::Ieee64}, data_value::DataValue};

enum Argument<'t> { Value(DataValue), Ref(&'t [u8]), Mut(&'t mut[u8]) }

use cranelift_interpreter::{interpreter::{InterpreterState, Interpreter}, step::Extension};
fn interpret(function: &Function, arguments: &mut [Argument], extension: Extension<DataValue>) {
	use Argument::*;
	let mut heap = vec![];
	let arguments_values = iter::map(arguments.iter(), |v| {
		match v {
			Value(v) => v.clone(),
			Ref(v) => { let base = heap.len(); heap.extend(v.iter()); DataValue::I32(base as i32) }
			Mut(v) => { let base = heap.len(); heap.extend(v.iter()); DataValue::I32(base as i32) }
		}
	});
	let mut interpreter = Interpreter::with_extension(InterpreterState{heap, ..default()}, extension);
	eprintln!("{:?}", &from_bytes::<f64>(&interpreter.state.heap));
	assert!(interpreter.call(function, &arguments_values).unwrap().unwrap_return() == &[]);
	eprintln!("{:?}", &from_bytes::<f64>(&interpreter.state.heap));
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

use cranelift_interpreter::value::ValueResult;
use combustion::reaction::Intrinsic;
pub fn extension(id: u32, arguments: &[DataValue]) -> ValueResult<DataValue> {
	let x = if let [DataValue::F64(x)] = arguments { x } else { unreachable!() };
	let x = f64::from_bits(x.bits());
	assert!(x.is_finite(), "{}", x);
	use std::convert::TryFrom;
	use Intrinsic::*;
	let op = Intrinsic::try_from(id).unwrap();
	let y = match op {
		exp2 => f64::exp2(x),
		log2 => f64::log2(x),
	};
	assert!(y.is_finite(), "{:?} {} = {}", op, x, y);
	Ok(DataValue::F64(Ieee64::/*from*/with_float(y)))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
	pretty_env_logger::init();
	let model = &std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let model = model::Model::new(&model)?;
	let (species_names, species) = combustion::Species::new(&model.species);
	use reaction::*;
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(&species_names, r));
	let width = 1;
	let states_len = ((1)/width)*width;
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
	let mut rates = vec![f64::NAN; (1/*2*/+species.len()-1)*states_len];

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
				Value(F64(Ieee64::with_float(pressure_R * reference_state.pressure_R))),
				Ref(as_bytes(temperature)),
				Ref(as_bytes(mass_fractions)),
				Value(F64(Ieee64::with_float(reference_state.temperature))),
				Mut(as_bytes_mut(mass_production_rates)),
				Value(F64(Ieee64::with_float(mass_production_rate_factor))),
				Value(F64(Ieee64::with_float(heat_release_rate_factor))),
				Mut(as_bytes_mut(heat_release_rate)),
			], extension);
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
    let rates = all_same(&rates, states_len);
		for ((&mass_production_rate, name), _molar_mass) in rates[1..].iter().zip(species_names.iter()).zip(species.molar_mass.iter()) {
			if mass_production_rate != 0. { println!("{:5} {:e}", name, mass_production_rate); }
		}
		println!("{:5} {:e}", "HRR", rates[0]);
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
