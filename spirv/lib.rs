#![feature(format_args_capture, default_free_fn)]
use {std::default::default, ast::*, itertools::Itertools, anyhow::Result};

struct FunctionBuilder {
	uniform:usize,
	input: usize,
	instructions: Vec<String>,
	values: Vec<Value>,
	results: usize,
}

impl FunctionBuilder {
fn expr(&mut self, e: &Expression) -> String {
	use Expression::*;
	match e {
		&F64(v) => {
			if v == f64::floor(v) {
				/*if v >= 1e4 { format!("{:e}", v) }
				else {*/ (v as i64).to_string()+"." //}
			}
			else { float_pretty_print::PrettyPrintFloat(v).to_string() }
		}
		Value(v) => {
			if v.0 < self.uniform{ format!("uniform[{}]", v.0) }
			else if v.0 < self.input { format!("input{}[id]", v.0) }
			else { assert!(self.values.contains(v), "{} {:?} {:?}", self.instructions.iter().format("\n"), self.values, v); format!("_{}", v.0) }
		},
		Neg(x) => format!("-{}", self.expr(x)),
		Max(a, b) => format!("max({}, {})", self.expr(a), self.expr(b)),
		Add(a, b) => format!("({} + {})", self.expr(a), self.expr(b)),
		Sub(a, b) => format!("({} - {})", self.expr(a), self.expr(b)),
		LessOrEqual(a, b) => format!("{} <= {}", self.expr(a), self.expr(b)),
		Mul(a, b) => format!("({} * {})", self.expr(a), self.expr(b)),
		Div(a, b) => format!("({} / {})", self.expr(a), self.expr(b)),
		Sqrt(x) => format!("sqrt({})", self.expr(x)),
		Block { statements, result } => {
			for s in &**statements { self.push(s) }
			let result = self.expr(result);
			let id = {let id = self.results; self.results += 1; id};
			self.instructions.push(format!("let _{id} = {result};"));
			format!("_{id}")
		}
		_ => panic!("{e:?}")
	}
}
fn push(&mut self, s: &Statement) {
	use Statement::*;
	match s {
		Value { id, value } => {
			let value = self.expr(value);
			self.instructions.push(format!("let _{} = {value};", id.0));
			assert!(!self.values.contains(id));
			self.values.push(id.clone());
		}
		//Output { index, value } => { let value = self.expr(value); self.instructions.push(format!("output{index}[id] = {value};")) }
		Branch { condition, consequent, alternative, results } => {
			self.instructions.extend(results.iter().map(|id| format!("var _{}: f64;", id.0)));
			self.values.extend(results.iter().cloned());
			let condition = self.expr(condition);
			self.instructions.push(format!("if ({condition}) {{"));
			for (id, value) in results.iter().zip(&**consequent) {
				let value = self.expr(value);
				self.instructions.push(format!("_{} = {value};", id.0));
			}
			self.instructions.push("} else {".into());
			for (id, value) in results.iter().zip(&**alternative) {
				let value = self.expr(value);
				self.instructions.push(format!("_{} = {value};", id.0));
			}
			self.instructions.push("}".into());
		}
	}
}
}

pub fn from(uniform: usize, Function{input, output, /*variables,*/ statements}: &Function) -> Result<Box<[u32]>> {
	let mut bindings = vec![format!("[[group(0), binding(0)]] var<push_constant> uniform : array<f64, {uniform}>;")];
	bindings.extend((uniform..*input).map(|i| { let binding = 1+i-uniform; format!("[[group(0), binding({binding})]] var<storage> input{i}: [[access(read)]] array<f64>;") }));
	bindings.extend((0..output.len()).map(|i| { let binding = 1+input-uniform+i; format!("[[group(0), binding({binding})]] var<storage> output{i}: [[access(write)]] array<f64>;") }));
	let module = [
		bindings,
		vec!["[[stage(compute), workgroup_size(1)]] fn main([[builtin(global_invocation_id)]] id3: vec3<u32>) { let id = id3.x;".to_string()],// var _: array<f64, {}>;", variables)],
		{
			let mut f = FunctionBuilder{uniform, input: *input, instructions: vec![], values: vec![], /*variables: vec![],*/ results: 0};
			for s in &**statements { f.push(s); }
			f.instructions
		},
		vec!["}".into()]
	].concat().join("\n");
	use naga::{front::wgsl::parse_str, valid::{Validator, Capabilities}, back::spv::{write_vec, Options, WriterFlags}};
	std::fs::write("/var/tmp/wgsl", &module)?;
	let module = parse_str(&module)?;
	Ok(write_vec(&module, &Validator::new(default(), Capabilities::all()).validate(&module)?, &Options{lang_version: (1, 6), flags: WriterFlags::DEBUG, ..default()})?.into())
	/*use {spirv::*, rspirv::dr::*}
	let mut b = Builder::new();
	b.set_version(1, 5);
	b.capability(Capability::Shader);
	b.memory_model(AddressingModel::Logical, MemoryModel::Vulkan);
	let return_type = b.type_void();
	let function_type = b.type_function(return_type, vec![]);
	let f = b.begin_function(return_type, /*id: */None, FunctionControl::DONT_INLINE | FunctionControl::CONST, function_type)?;
	b.begin_basic_block(None).unwrap();
	b.ret().unwrap();
	b.end_function().unwrap();
	b.entry_point(spirv::ExecutionModel::GLCompute, f, "main", vec![]);
	b.module().assemble().into()*/
}
