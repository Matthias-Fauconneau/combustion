#![feature(default_free_fn, in_band_lifetimes, bindings_after_at, unboxed_closures)] #![allow(non_snake_case)]
use {fehler::throws, anyhow::Error};

pub fn as_bytes<T>(value: &T) -> &[u8] { unsafe{std::slice::from_raw_parts(value as *const T as *const u8, std::mem::size_of::<T>())} }

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

mod vulkan;

use cranelift_codegen::ir::{Function, AbiParam, types::{F64, I32}};

#[throws(std::fmt::Error)] fn glsl(function: Function) -> String {
	use std::convert::TryFrom;
	let mut w = String::from(
	r#"#version 460
double fneg(double x) { return -x; }
double fmax(double x, double y) { return max(x, y); }
double fadd(double x, double y) { return x * y; }
double fsub(double x, double y) { return x - y; }
double fmul(double x, double y) { return x * y; }
double fdiv(double x, double y) { return x / y; }
uint ishl_imm(uint x, uint y) { return x << y; }
uint iadd(uint x, uint y) { return x + y; }
layout(set = 0, binding = 0) readonly restrict buffer Load { double[] load; } load;
layout(set = 0, binding = 1) restrict buffer Store { double[] store; } store;
layout(push_constant) uniform Constants {
"#);
	use std::fmt::Write;
	use itertools::Itertools;
	let to_string = |value_type| {
		match value_type {
			F64 => "double",
			I32 => "uint",
			_ => unimplemented!(),
		}
	};
	write!(w, "{}", function.signature.params.iter().enumerate().skip(1).map(|(i, AbiParam{value_type, ..})| format!("{} v{};\n", to_string(*value_type), i)).format(""))?;
	w.push_str(r#"} constants; void main() {
	const uint v0 = gl_GlobalInvocationID.x;
"#);
	write!(w, "{}", function.signature.params.iter().enumerate().skip(1).map(
		|(i, AbiParam{value_type, ..})| format!("{} v{i} = constants.v{i};\n", to_string(*value_type), i=i)
	).format(""))?;
	let f = function.dfg;
	use iter::Single;
	for instruction in function.layout.block_insts(function.layout.blocks().single().unwrap()) {
		match f.inst_results(instruction) {
			[] => (),
			[result] => write!(w, "{} {} = ", to_string(f.value_type(*result)), result)?,
			_ => unimplemented!(),
		};
		use cranelift_codegen::ir::{InstructionData::*, instructions::Opcode};
		use combustion::reaction::Intrinsic;
		fn i32_from(offset: impl Into<i32>) -> i32 { offset.into() } // Workaround impl Into !From
		match &f[instruction] {
			UnaryIeee64{imm, ..} => write!(w, "{}", {let s = f64::from_bits(imm.bits()).to_string(); if s.contains('.') { s } else { s+"." }})?,
			Unary{opcode, arg} => write!(w, "{}({})", opcode, arg)?,
			Load{arg, offset, ..} => write!(w, "load.load[({}+{})/8]", arg, i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/)?,
			Store{args, offset, ..} => write!(w, "store.store[({}+{})/8] = {}", args[1], i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/, args[0])?,
			Call{func_ref, args, ..} => write!(w, "{:?}({})", Intrinsic::try_from(func_ref.as_u32()).unwrap(), args.first(&f.value_lists).unwrap())?,
			Binary{opcode, args} => write!(w, "{}({}, {})", opcode, args[0], args[1])?,
			BinaryImm64{opcode, arg, imm} => write!(w, "{}({}, {})", opcode, arg, imm)?,
			MultiAry{opcode: Opcode::Return, ..} => write!(w, "return")?,
			instruction => unimplemented!("{:?}", instruction)
		};
		write!(w, ";\n")?;
	}
	write!(w, "}}")?;
	w
}

use reaction::Simulation;
#[throws] fn main() {
	let model = std::fs::read("CH4+O2.ron")?;
	let simulation = Simulation::new(&model)?;
	let states_len = simulation.states_len();
	let Simulation{species_names, function, states, rates, pressure_Pa_R, reference_temperature, mass_production_rate_factor, heat_release_rate_factor} = simulation;
	let function = glsl(function)?;
	std::fs::write("/var/tmp/main.comp", &function)?;
  let main = shaderc::Compiler::new().unwrap().compile_into_spirv(&function, shaderc::ShaderKind::Compute, "main.comp", "main", None)?;

	use vulkan::{Device, Buffer};
	let ref device = Device::new()?;
	let load = Buffer::new(device, states.iter().copied())?;
	let store = Buffer::new(device, rates.iter().copied())?;
	let constants = reaction::Constants{
		pressure_Pa_R,
		temperature: 0,
		mass_fractions: (states_len*std::mem::size_of::<f64>()) as u32,
		mass_production_rates: (states_len*std::mem::size_of::<f64>()) as u32,
		heat_release_rate: 0,
		reference_temperature,
		mass_production_rate_factor,
		heat_release_rate_factor
	};
	let buffers = [&[&load] as &[_], &[&store] as &[_]];
	let pipeline = time!(device.pipeline(main.as_binary(), as_bytes(&constants), &buffers))?; // Compiles SPIRV -> Gen
	device.bind(pipeline.descriptor_set, &buffers)?;
	let width = 1;
	let command_buffer = device.command_buffer(&pipeline, as_bytes(&constants), width, states_len)?;
	dbg!();
	let time = device.submit_and_wait(command_buffer)?;
	dbg!();
	println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e3, (states_len as f32)/1e3/time);
	let rates = store.map(device).unwrap();
	reaction::report(&species_names, &rates);
}
