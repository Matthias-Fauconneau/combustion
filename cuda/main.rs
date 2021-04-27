#![feature(bool_to_option)]#![allow(non_snake_case)]
use {fehler::throws, anyhow::Error};
use cranelift_codegen::ir::{Function, AbiParam, types::{F64, I32}};

#[throws(std::fmt::Error)] fn cu(function: Function) -> String {
	use std::convert::TryFrom;
	let mut w = String::from(
	r#"__device__ double neg(double x) { return -x; }
	__device__ double add(double x, double y) { return x*y; }
	__device__ double sub(double x, double y) { return x-y; }
	__device__ double mul(double x, double y) { return x*y; }
	__device__ double div(double x, double y) { return x/y; }
	__global__ void kernel(
	"#);
	use std::fmt::Write;
	use itertools::Itertools;
	write!(w, "{}", function.signature.params.iter().enumerate().skip(1).map(|(i, AbiParam{value_type, ..})|
		match *value_type {
			F64 => "double",
			I32 => "double*",
			_ => unimplemented!(),
		}.to_owned() + &format!(" p{}", i)
	).format(", "))?;
	w.push_str(r#") {
	const uint v0 = blockIdx.x * /*SIMD width*/blockDim.x + /*SIMD lane*/threadIdx.x;
"#);
	for (i, AbiParam{value_type, ..}) in function.signature.params.iter().enumerate().skip(1) {
		match *value_type {
			F64 => write!(w, "double v{i} = p{i};\n", i=i)?,
			I32 => write!(w, "double* v{i} = p{i}+v0;\n", i=i)?,
			_ => unimplemented!(),
		};
	}
	let f = function.dfg;
	use iter::Single;
	for instruction in function.layout.block_insts(function.layout.blocks().single().unwrap()) {
		match f.inst_results(instruction) {
			[] => (),
			[result] => write!(w, "double {} = ", result)?,
			_ => unimplemented!(),
		};
		use cranelift_codegen::ir::InstructionData::*;
		use combustion::reaction::Intrinsic;
		fn i32_from(offset: impl Into<i32>) -> i32 { offset.into() } // Workaround impl Into !From
		match &f[instruction] {
			UnaryIeee64{imm, ..} => write!(w, "{}", {let s = f64::from_bits(imm.bits()).to_string(); if s.contains('.') { s } else { s+"." }})?,
			Unary{opcode, arg} => write!(w, "{}({})", opcode, arg)?,
			Load{arg, offset, ..} => write!(w, "{}[{}]", arg, i32_from(*offset)/(std::mem::size_of::<f64>() as i32))?,
			Store{args, offset, ..} => write!(w, "{}[{}] = {}", args[1], i32_from(*offset)/(std::mem::size_of::<f64>() as i32), args[0])?,
			Call{func_ref, args, ..} => write!(w, "{:?}({})", Intrinsic::try_from(func_ref.as_u32()).unwrap(), args.first(&f.value_lists).unwrap())?,
			Binary{opcode, args} => write!(w, "{}({}, {})", opcode, args[0], args[1])?,
			instruction => unimplemented!("{:?}", instruction)
		};
		write!(w, ";\n")?;
	}
	write!(w, "}}")?;
	w.replace("fneg", "neg").replace("fadd", "add").replace("fsub", "sub").replace("fmul", "mul").replace("fdiv", "div")
}

#[throws] fn main() {
	let model = std::fs::read("CH4+O2.ron")?;
	use reaction::Simulation;
	let simulation = Simulation::new(&model)?;
	let states_len = simulation.states_len();
	let Simulation{species_names, function, states, rates, pressure_Pa_R, reference_temperature, mass_production_rate_factor, heat_release_rate_factor} = simulation;
	std::fs::write("/var/tmp/main.cu", &cu(function)?)?;
	std::process::Command::new("nvcc").args(&["--ptx","/var/tmp/main.cu","-o","/var/tmp/main.ptx"]).spawn()?.wait()?.success().then_some(()).unwrap();
	rustacuda::init(CudaFlags::empty()).unwrap();
	use rustacuda::{prelude::*, launch};
	let device = Device::get_device(0).unwrap();
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).unwrap();
	let module = Module::load_from_string(&std::ffi::CString::new(std::fs::read("/var/tmp/main.ptx").unwrap()).unwrap()).unwrap();
	let stream = Stream::new(StreamFlags::NON_BLOCKING, None).unwrap();
	let mut device_states = DeviceBuffer::from_slice(&states).unwrap();
	let mut device_rates = DeviceBuffer::from_slice(&rates).unwrap();
	let width = 1;
	let start = std::time::Instant::now();
	unsafe{launch!(module._Z6kerneldPdS_S_ddS_<<</*workgroupCount*/(states_len/width) as u32,/*workgroupSize*/width as u32, 0, stream>>>(
		pressure_Pa_R,
		device_states.as_device_ptr(),
		reference_temperature,
		device_states.as_device_ptr().add(states_len),
		device_rates.as_device_ptr().add(states_len),
		mass_production_rate_factor,
		heat_release_rate_factor,
		device_rates.as_device_ptr()
	)).unwrap()}
	stream.synchronize().unwrap();
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f64();
	println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e3, (states_len as f64)/1e3/time);
	let mut rates = rates;
	device_rates.copy_to(&mut rates).unwrap();
	reaction::report(&species_names, &rates);
}
