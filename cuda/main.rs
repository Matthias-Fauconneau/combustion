#![feature(bool_to_option)]#![allow(non_snake_case)]
use {fehler::throws, anyhow::Error};
use cranelift_codegen::ir::{Function, AbiParam, types::{F64, I32}};

fn test() -> Function {
	let mut function = Function::new();
	let mut function_builder_context = cranelift_frontend::FunctionBuilderContext::new();
	let mut f = cranelift_frontend::FunctionBuilder::new(&mut function, &mut function_builder_context);
	let entry_block = f.create_block();
	let parameters = [I32, F64, I32, I32, I32, I32, F64, F64, F64];
	f.func.signature.params = iter::map(&parameters, |t| AbiParam::new(*t)).to_vec();
	f.append_block_params_for_function_params(entry_block);
	f.switch_to_block(entry_block);
	f.seal_block(entry_block);
	use std::convert::TryInto;
	use cranelift_codegen::ir::{InstBuilder, entities::Value, MemFlags, immediates::Ieee64};
	let [index, _constant, /*state*/T, mass_fractions, /*rates*/mass_production_rates, heat_release_rates,
				_reference_temperature, _mass_production_rates_factor, _heat_release_rate_factor]: [Value; 9] = f.block_params(entry_block).try_into().unwrap();
	let offset = f.ins().ishl_imm(index, 3);
	let T = f.ins().iadd(T, offset);
	let _mass_fractions = f.ins().iadd(mass_fractions, offset);
	let mass_production_rates = f.ins().iadd(mass_production_rates, offset);
	let heat_release_rates = f.ins().iadd(heat_release_rates, offset);
	let T = f.ins().load(F64, MemFlags::trusted(), T, 0);
	f.ins().store(MemFlags::trusted(), T, heat_release_rates, 0);
	for i in 0..52 {
		let c = f.ins().f64const(Ieee64::with_float(i as f64));
		f.ins().store(MemFlags::trusted(), c, mass_production_rates, (i*std::mem::size_of::<f64>() as u64) as i32);
	}
	f.ins().return_(&[]);
	function
}

#[throws(std::fmt::Error)] fn cu(function: Function) -> String {
	use std::convert::TryFrom;
	let mut w = String::from(
	r#"__device__ double neg(double x) { return -x; }
	__device__ double add(double x, double y) { return x * y; }
	__device__ double sub(double x, double y) { return x - y; }
	__device__ double mul(double x, double y) { return x * y; }
	__device__ double div(double x, double y) { return x / y; }
	__device__ unsigned long ishl_imm(unsigned long x, unsigned int y) { return x << y; }
	__device__ unsigned long iadd(unsigned long x, unsigned int y) { return x + y; }
	__global__ void kernel(
	"#);
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
	w.push_str(r#") {
	const unsigned int v0 = blockIdx.x * /*SIMD width*/blockDim.x + /*SIMD lane*/threadIdx.x;
"#);
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
			Load{arg, offset, ..} => write!(w, "*(double*)({}+{})", arg, i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/)?,
			Store{args, offset, ..} => write!(w, "*(double*)({}+{}) = {}", args[1], i32_from(*offset)/*(std::mem::size_of::<f64>() as i32)*/, args[0])?,
			Call{func_ref, args, ..} => write!(w, "{:?}({})", Intrinsic::try_from(func_ref.as_u32()).unwrap(), args.first(&f.value_lists).unwrap())?,
			Binary{opcode, args} => write!(w, "{}({}, {})", opcode, args[0], args[1])?,
			BinaryImm64{opcode, arg, imm} => write!(w, "{}({}, {})", opcode, arg, imm)?,
			MultiAry{opcode: Opcode::Return, ..} => write!(w, "return")?,
			instruction => unimplemented!("{:?}", instruction)
		};
		write!(w, ";\n")?;
	}
	write!(w, "}}")?;
	w.replace("fneg", "neg").replace("fadd", "add").replace("fsub", "sub").replace("fmul", "mul").replace("fdiv", "div")
}

use reaction::Simulation;
#[throws] fn main() {
	let model = std::fs::read("CH4+O2.ron")?;
	let simulation = Simulation::new(&model)?;
	let states_len = simulation.states_len();
	let Simulation{species_names, function, states, rates, pressure_Pa_R, reference_temperature, mass_production_rate_factor, heat_release_rate_factor} = simulation;
	std::fs::write("/var/tmp/main.cu", &cu(function)?)?;
	//std::fs::write("/var/tmp/main.cu", &cu(test())?)?;
	std::process::Command::new("nvcc").args(&["--ptx","/var/tmp/main.cu","-o","/var/tmp/main.ptx"]).spawn()?.wait()?.success().then_some(()).unwrap();
	rustacuda::init(CudaFlags::empty()).unwrap();
	use rustacuda::{prelude::*, launch};
	let device = Device::get_device(0).unwrap();
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device).unwrap();
	let module = Module::load_from_string(&std::ffi::CString::new(std::fs::read("/var/tmp/main.ptx").unwrap()).unwrap()).unwrap();
	let stream = Stream::new(/*irrelevant*/StreamFlags::NON_BLOCKING, None).unwrap();
	let mut device_states = DeviceBuffer::from_slice(&states).unwrap();
	let mut device_rates = DeviceBuffer::from_slice(&rates).unwrap();
	let width = 1;
	let start = std::time::Instant::now();
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
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f64();
	println!("{:.0}K in {:.1}ms = {:.2}ms, {:.1}K/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e3, (states_len as f64)/1e3/time);
	let mut rates = rates;
	device_rates.copy_to(&mut rates).unwrap();
	reaction::report(&species_names, &rates);
}
