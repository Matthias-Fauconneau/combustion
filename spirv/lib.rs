#![feature(format_args_capture, default_free_fn)]
use std::default::default;
use ast::*;
use itertools::Itertools;

struct WGSL {
	instructions: Vec<String>,
	definitions: Vec<Value>,
	variables: Vec<Variable>,
	results: usize,
}

impl WGSL {
fn expr(&mut self, e: &Expression) -> String {
	use Expression::*;
	match e {
		&Literal(v) => {
			if v == f64::floor(v) {
				if v >= 1e4 { format!("{:e}", v) }
				else { (v as i64).to_string() }
			}
			else { float_pretty_print::PrettyPrintFloat(v).to_string() }
		}
		Use(v) => { assert!(self.definitions.contains(v)); format!("_{}", v.0) },
		Load(v) => { assert!(self.variables.contains(v)); format!("__[{}]", v.0) },
		Index { base, index } => format!("arrays[{}][{index}]", base.0),
		Neg(x) => format!("-{}", self.expr(x)),
		Add(a, b) => format!("({} + {})", self.expr(a), self.expr(b)),
		Sub(a, b) => format!("({} - {})", self.expr(a), self.expr(b)),
		Less(a, b) => format!("{} < {}", self.expr(a), self.expr(b)),
		Mul(a, b) => format!("({} * {})", self.expr(a), self.expr(b)),
		Div(a, b) => format!("({} / {})", self.expr(a), self.expr(b)),
		Call { function, arguments } => format!("{function}({})", arguments.iter().map(|e| self.expr(e)).format(", ")),
		Block { statements, result } => {
			for s in &**statements { self.push(s) }
			let result = self.expr(result);
			let id = {let id = self.results; self.results += 1; id};
			self.instructions.push(format!("let _{id} = {result};"));
			"_{id}".into()
		}
	}
}
fn push(&mut self, s: &Statement) {
	use Statement::*;
	match s {
		Define { id, value } => {
			let value = self.expr(value);
			self.instructions.push(format!("let _{} = {value};", id.0));
			assert!(!self.definitions.contains(id));
			self.definitions.push(Value(id.0));
		}
		Store { id, value } => { let value = self.expr(value); self.instructions.push(format!("__[{}] = {value};", id.0)) }
		Output { index, value } => { let value = self.expr(value); self.instructions.push(format!("_[{index}] = {value};")) }
		Branch { condition, consequent, alternative } => {
			let condition = self.expr(condition);
			self.instructions.push(format!("if ({condition}) {{"));
			for s in &**consequent { self.push(s) }
			self.instructions.push("} else {".into());
			for s in &**alternative { self.push(s) }
		}
	}
}
}

pub fn from<const U: usize, const V: usize, const A: usize>(function: &Function<U,V,A>) -> Result<Box<[u32]>, Box<dyn std::error::Error>> {
	use naga::{front::wgsl::parse_str, valid::{Validator, Capabilities}, back::spv::{write_vec, Options, WriterFlags}};
	let module = parse_str(&format!(r#"
		[[set(0), binding(0)]] var<push_constant> parameters : array<f64, {V}>;
		[[set(1), binding(1)]] var<storage,readonly> arrays: array<array<f64>>
		[[set(2), binding(2)]] var<storage,writeonly> _: array<f64>
		[[stage(compute), workgroup_size(1)]]
		fn main([[builtin(global_invocation_id)]] id: vec3<u32>) {{
			var _[{}];
			var __[{}];
			{}
			return _;
		}}"#, function.output, function.variables, {
			let mut wgsl = WGSL{instructions: vec![], definitions: vec![], variables: vec![], results: 0};
			for s in &*function.statements { wgsl.push(s); }
			wgsl.instructions.join("\n")
		}
	))?;
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
