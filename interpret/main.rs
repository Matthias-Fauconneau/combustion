#![allow(non_snake_case)]#![feature(bool_to_option,assert_matches,default_free_fn)]
use {std::default::default, fehler::throws, anyhow::Error};

pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
pub fn as_bytes_mut<T>(slice: &mut [T]) -> &mut [u8] { unsafe{std::slice::from_raw_parts_mut(slice.as_mut_ptr() as *mut u8, slice.len() * std::mem::size_of::<T>())} }
pub fn from_bytes<T>(slice: &[u8]) -> &[T] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const T, slice.len() / std::mem::size_of::<T>())} }

use cranelift_codegen::{ir::{function::Function, immediates::Ieee64, Value}, data_value::DataValue};

enum Argument<'t> { Value(DataValue), Ref(&'t [u8]), Mut(&'t mut[u8]) }

use cranelift_interpreter::{interpreter::{InterpreterState, Interpreter}, step::Extension};
fn interpret(function: &Function, arguments: &mut [Argument], extension: Extension<DataValue>) -> Box<[Option<DataValue>]> {
	let mut heap = vec![];
	let arguments_values = iter::map(arguments.iter(), |v| {
		use Argument::*;
		match v {
			Value(v) => v.clone(),
			Ref(v) => { let base = heap.len(); heap.extend(v.iter()); DataValue::I32(base as i32) }
			Mut(v) => { let base = heap.len(); heap.extend(v.iter()); DataValue::I32(base as i32) }
		}
	});
	let mut interpreter = Interpreter::with_extension(InterpreterState{heap, ..default()}, extension);
	let (_, values) = interpreter.call(function, &arguments_values).unwrap().unwrap_return();
	let values = values.into_boxed_slice();
	let mut base = 0;
	let state = interpreter.state;
	for v in arguments.iter_mut() {
		use Argument::*;
		match v {
			Value(_) => (),
			Ref(v) => base += v.len(),
			Mut(v) => { v.copy_from_slice(&state.heap[base..base+v.len()]); base += v.len(); }
		}
	}
	values
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

use reaction::Simulation;
#[throws] fn main() {
	pretty_env_logger::init();
	let model = std::fs::read("CH4+O2.ron")?;
	let simulation = Simulation::new(&model)?;
	let states_len = simulation.states_len();
	let Simulation{species_names, function, states, mut rates, pressure_Pa_R, reference_temperature, mass_production_rate_factor, heat_release_rate_factor} = simulation;
	let (temperature, mass_fractions) = states.split_at(states_len);
	let (heat_release_rate, mass_production_rates) = rates.split_at_mut(states_len);
	use DataValue::*;
	let start = std::time::Instant::now();
	let values = interpret(&function, {use Argument::*; &mut [
		Value(I32(0)),
		Value(F64(Ieee64::with_float(pressure_Pa_R))),
		Ref(as_bytes(temperature)),
		Ref(as_bytes(mass_fractions)),
		Mut(as_bytes_mut(mass_production_rates)),
		Mut(as_bytes_mut(heat_release_rate)),
		Value(F64(Ieee64::with_float(reference_temperature))),
		Value(F64(Ieee64::with_float(mass_production_rate_factor))),
		Value(F64(Ieee64::with_float(heat_release_rate_factor))),
	]}, extension);
	let other_values = iter::map(from_bytes(&std::fs::read("/var/tmp/values")?), |value:&f64| if value.is_nan() { None } else { Some(DataValue::F64(Ieee64::with_float(*value))) });
	/*for (ref_value, value) in ref_values.iter().zip(values.iter()) {
		println!("{:?} {:?}", ref_value, value);
		match (ref_value, value) {
			(Some(I32(_)), None) => {},
			(Some(F64(_)), None) => {},
			_ => assert_eq!(ref_value, value),
		}
	}*/
	let value = |value:Value| {
		let value = value.as_u32() as usize;
		fn to_string(x: &DataValue) -> String {
			match x {
				F64(x) => format!("{}", f64::from_bits(x.bits())),
				x => format!("{}", x),
			}
		}
		match (&values[value], &other_values[value]) {
			(Some(v), None) => format!("{}", to_string(v)),
			(Some(a), Some(b)) => if a == b { format!("{}", to_string(a)) } else { format!("{}/{}", to_string(a), to_string(b)) },
			_ => unimplemented!(),
		}
	};
	let f = function.dfg;
	use iter::Single;
	for instruction in function.layout.block_insts(function.layout.blocks().single().unwrap()) {
		match f.inst_results(instruction) {
			[] => (),
			[result] => print!("#{:>3} = ", result.as_u32()),
			_ => unimplemented!(),
		};
		use cranelift_codegen::ir::{InstructionData::*, instructions::Opcode};
		use combustion::reaction::Intrinsic;
		use std::convert::TryFrom;
		fn i32_from(offset: impl Into<i32>) -> i32 { offset.into() } // Workaround impl Into !From
		print!("{:>60}", match &f[instruction] {
			UnaryIeee64{imm, ..} => format!("{}", {let s = f64::from_bits(imm.bits()).to_string(); if s.contains('.') { s } else { s+"." }}),
			Unary{opcode, arg} => format!("{}({}={})", opcode, arg, value(*arg)),
			Load{arg, offset, ..} => format!("[{}+{}]", arg, i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/),
			Store{args, offset, ..} => format!("[{}+{}] = {}", args[1], i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/, args[0]),
			Call{func_ref, args, ..} => format!("{:?}({})", Intrinsic::try_from(func_ref.as_u32()).unwrap(), args.first(&f.value_lists).unwrap()),
			Binary{opcode, args} => format!("{}({}={}, {}={})", opcode, args[0], value(args[0]), args[1], value(args[1])),
			BinaryImm64{opcode, arg, imm} => format!("{}({}, {})", opcode, arg, imm),
			MultiAry{opcode: Opcode::Return, ..} => format!("return"),
			instruction => unimplemented!("{:?}", instruction)
		});
		match &f[instruction] {
			UnaryIeee64{..} => {},
			_ => match f.inst_results(instruction) {
				[] => (),
				[result] => {
					print!(" => {}", value(*result));
					let result = result.as_u32() as usize;
					match (&values[result], &other_values[result]) {
						(Some(_), None) => {},
						(a, b) => assert_eq!(a, b),
					}
				}
				_ => unimplemented!(),
			},
		};
		print!("\n");
	}
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f64();
	println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e3, (states_len as f64)/1e3/time);
	reaction::report(&species_names, &rates);
}
