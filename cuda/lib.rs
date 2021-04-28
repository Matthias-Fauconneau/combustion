use cranelift_codegen::ir::{Function, AbiParam, types::{F64, I32}};

#[fehler::throws(std::fmt::Error)] pub fn cu(function: Function, store_values: bool) -> String {
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
