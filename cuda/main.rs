#![feature(bool_to_option)]#![allow(non_snake_case)]
use {fehler::throws, anyhow::Error};
use cranelift_codegen::ir::{Function, AbiParam, types::{F64, I32}};

#[throws(std::fmt::Error)] fn cu(function: Function, store_values: bool) -> String {
	use std::convert::TryFrom;
	let mut w = String::from("__global__ void kernel(");
	w.reserve(1024*1024);
	use std::fmt::Write;
	use itertools::Itertools;
	let to_string = |value_type| {
		match value_type {
			F64 => "double",
			I32 => "unsigned long",
			_ => unimplemented!(),
		}
	};
	write!(w, "{}", function.signature.params.iter().enumerate().skip(1).map(|(i, AbiParam{value_type, ..})| format!("{} v{}", to_string(*value_type), i)).format(", "))?;
	if store_values { w.push_str(", double* values"); }
	w.push_str(r#") {
	const unsigned int v0 = blockIdx.x * /*SIMD width*/blockDim.x + /*SIMD lane*/threadIdx.x;
"#);
	if store_values { for (i, AbiParam{value_type, ..}) in function.signature.params.iter().enumerate().skip(1) { if *value_type == F64 { write!(w, "values[{i}] = v{i};\n", i=i)? } } }
	let f = function.dfg;
	use iter::Single;
	for instruction in function.layout.block_insts(function.layout.blocks().single().unwrap()) {
		match f.inst_results(instruction) {
			[] => (),
			[result] => write!(w, "{} {} = ", to_string(f.value_type(*result)), result)?,
			_ => unimplemented!(),
		};
		use cranelift_codegen::ir::{InstructionData::{self, *}, instructions::Opcode::*};
		use combustion::reaction::Intrinsic;
		fn i32_from(offset: impl Into<i32>) -> i32 { offset.into() } // Workaround impl Into !From
		match &f[instruction] {
			UnaryIeee64{imm, ..} => write!(w, "{}", {let s = f64::from_bits(imm.bits()).to_string(); if s.contains('.') { s } else { s+"." }})?,
			//Unary{opcode, arg} => write!(w, "{}({})", opcode, arg)?,
			Unary{opcode: Fneg, arg} => write!(w, "-{}", arg)?,
			InstructionData::Load{arg, offset, ..} => write!(w, "*(double*)({}+{})", arg, i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/)?,
			InstructionData::Store{args, offset, ..} => write!(w, "*(double*)({}+{}) = {}", args[1], i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/, args[0])?,
			InstructionData::Call{func_ref, args, ..} => write!(w, "{:?}({})", Intrinsic::try_from(func_ref.as_u32()).unwrap(), args.first(&f.value_lists).unwrap())?,
			//Binary{opcode, args} => write!(w, "{}({}, {})", opcode, args[0], args[1])?,
			Binary{opcode: Iadd, args} => write!(w, "{}+{}", args[0], args[1])?,
			Binary{opcode: Fmax, args} => write!(w, "fmax({}, {})", args[0], args[1])?,
			Binary{opcode: Fadd, args} => write!(w, "{}+{}", args[0], args[1])?,
			Binary{opcode: Fsub, args} => write!(w, "{}-{}", args[0], args[1])?,
			Binary{opcode: Fmul, args} => write!(w, "{}*{}", args[0], args[1])?,
			Binary{opcode: Fdiv, args} => write!(w, "{}/{}", args[0], args[1])?,
			//BinaryImm64{opcode, arg, imm} => write!(w, "{}({}, {})", opcode, arg, imm)?,
			BinaryImm64{opcode: IshlImm, arg, imm} => write!(w, "{} << {}", arg, imm)?,
			MultiAry{opcode: Return, ..} => write!(w, "return")?,
			instruction => unimplemented!("{:?}", instruction)
		};
		write!(w, ";\n")?;
		if store_values {
			if let [result ] = f.inst_results(instruction) { if f.value_type(*result) == F64 {
				assert!(result.as_u32() < f.values().count() as u32);
				write!(w, "values[{}] = {};", result.as_u32(), result)?
			} }
		}
	}
	write!(w, "}}")?;
	w
}

use reaction::Simulation;
#[throws] fn main() {
	let model = std::fs::read("CH4+O2.ron")?;
	let width = 32;
	let simulation = Simulation::new(&model, 512*width)?;
	let states_len = simulation.states_len();
	let Simulation{species_names, function, states, rates, pressure_Pa_R, reference_temperature, mass_production_rate_factor, heat_release_rate_factor} = simulation;
	let values_count = function.dfg.values().count();
	let store_values = false;
	let function = cu(function, store_values)?;
	if std::fs::metadata("/var/tmp/main.cu").map_or(true, |cu|
		std::fs::read("/var/tmp/main.cu").unwrap() != function.as_bytes() ||
		std::fs::metadata("/var/tmp/main.ptx").map_or(true, |ptx| ptx.modified().unwrap() < cu.modified().unwrap())
	) {
		std::fs::write("/var/tmp/main.cu", &function)?;
		std::process::Command::new("nvcc").args(&["--ptx","/var/tmp/main.cu","-o","/var/tmp/main.ptx"]).spawn()?.wait()?.success().then_some(()).unwrap();
	}
	rustacuda::init(CudaFlags::empty()).unwrap();
	use rustacuda::{prelude::*, launch};
	let device = Device::get_device(0).unwrap();
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).unwrap();
	let module = Module::load_from_string(&std::ffi::CString::new(std::fs::read("/var/tmp/main.ptx").unwrap()).unwrap()).unwrap();
	let stream = Stream::new(/*irrelevant*/StreamFlags::NON_BLOCKING, None).unwrap();
	let mut device_states = DeviceBuffer::from_slice(&states).unwrap();
	let mut device_rates = DeviceBuffer::from_slice(&rates).unwrap();
	let start = std::time::Instant::now();
	if store_values {
		let values = vec![f64::NAN; values_count];
		let mut device_values = DeviceBuffer::from_slice(&values).unwrap();
		unsafe{launch!(module._Z6kerneldmmmmdddPd<<</*workgroupCount*/(states_len/width) as u32,/*workgroupSize*/width as u32, 0, stream>>>(
			pressure_Pa_R,
			device_states.as_device_ptr(),
			device_states.as_device_ptr().add(states_len),
			device_rates.as_device_ptr().add(states_len),
			device_rates.as_device_ptr(),
			reference_temperature,
			mass_production_rate_factor,
			heat_release_rate_factor,
			device_values.as_device_ptr()
		)).unwrap()}
		stream.synchronize().unwrap();
		let mut values = values;
		device_values.copy_to(&mut values).unwrap();
		pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		std::fs::write("/var/tmp/values", as_bytes(&values))?;
	} else {
		unsafe{launch!(module._Z6kerneldmmmmddd<<</*workgroupCount*/(states_len/width) as u32,/*workgroupSize*/width as u32, 0, stream>>>(
			pressure_Pa_R,
			device_states.as_device_ptr(),
			device_states.as_device_ptr().add(states_len),
			device_rates.as_device_ptr().add(states_len),
			device_rates.as_device_ptr(),
			reference_temperature,
			mass_production_rate_factor,
			heat_release_rate_factor
		)).unwrap()}
		stream.synchronize().unwrap();
	}
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f64();
	println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e3, (states_len as f64)/1e3/time);
	let mut rates = rates;
	device_rates.copy_to(&mut rates).unwrap();
	reaction::report(&species_names, &rates);
}
