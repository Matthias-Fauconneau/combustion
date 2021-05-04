#![allow(non_snake_case)]#![feature(bool_to_option,assert_matches,default_free_fn)]
use {std::default::default, fehler::throws, anyhow::Error};

pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
pub fn as_bytes_mut<T>(slice: &mut [T]) -> &mut [u8] { unsafe{std::slice::from_raw_parts_mut(slice.as_mut_ptr() as *mut u8, slice.len() * std::mem::size_of::<T>())} }
pub fn from_bytes<T>(slice: &[u8]) -> &[T] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const T, slice.len() / std::mem::size_of::<T>())} }

use cranelift_codegen::{ir::{function::Function, immediates::Ieee64, Value}, data_value::DataValue};

enum Argument<'t> { Value(DataValue), Ref(&'t [u8]), Mut(&'t mut[u8]) }

use cranelift_interpreter::{interpreter::{InterpreterState, Interpreter}};
fn interpret(function: &Function, arguments: &mut [Argument]) -> Box<[Option<DataValue>]> {
	let mut heap = vec![];
	let arguments_values = iter::map(arguments.iter(), |v| {
		use Argument::*;
		match v {
			Value(v) => v.clone(),
			Ref(v) => { let base = heap.len(); heap.extend(v.iter()); DataValue::I32(base as i32) }
			Mut(v) => { let base = heap.len(); heap.extend(v.iter()); DataValue::I32(base as i32) }
		}
	});
	let mut interpreter = Interpreter::new(InterpreterState{heap, ..default()});
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

use reaction::Simulation;
#[throws] fn main() {
	pretty_env_logger::init();
	let model = std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let model = model::Model::new(&model)?;
	let (species_names, species) = combustion::Species::new(&model.species);
	use reaction::*;
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(&species_names, r));
	let function = rate::<_,{Property::Pressure}>(&species, &*reactions, 1);
	let ref state = initial_state(&model);
	let rates = vec![f64::NAN; (2+species.len()-1)*1].into_boxed_slice();
	let start = std::time::Instant::now();
	let values = interpret(&function, {use {DataValue::*, Argument::*}; &mut [
		Value(I32(0)),
		Value(F64(Ieee64::with_float(pressure_Pa_R))),
		Ref(as_bytes(states)),
		Mut(as_bytes_mut(rates)),
	]});
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
		use std::convert::TryFrom;
		fn i32_from(offset: impl Into<i32>) -> i32 { offset.into() } // Workaround impl Into !From
		print!("{:>60}", match &f[instruction] {
			UnaryIeee64{imm, ..} => format!("{}", {let s = f64::from_bits(imm.bits()).to_string(); if s.contains('.') { s } else { s+"." }}),
			Unary{opcode, arg} => format!("{}({}={})", opcode, arg, value(*arg)),
			Load{arg, offset, ..} => format!("[{}+{}]", arg, i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/),
			Store{args, offset, ..} => format!("[{}+{}] = {}", args[1], i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/, args[0]),
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
	reaction::report(&species_names, &rates);
}
