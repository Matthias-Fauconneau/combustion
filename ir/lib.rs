#![feature(default_free_fn,format_args_capture,in_band_lifetimes,array_map)]
use std::default::default;
pub use cranelift::codegen::ir::{function::Function, types::{I32, I64, F32}};
use cranelift::{
	frontend::{self, FunctionBuilderContext},
	codegen::{ir::{InstBuilder, AbiParam, entities::Value, MemFlags}, }, //types::Type,
};

#[derive(derive_more::Deref,derive_more::DerefMut)] struct FunctionBuilder<'t> {
	#[deref]#[deref_mut] builder: frontend::FunctionBuilder<'t>,
	constants_u32: std::collections::HashMap<u32, Value>,
	constants_f32: std::collections::HashMap<u32, Value>,
}

impl FunctionBuilder<'t> {
	fn new(function: &'t mut Function, function_builder_context: &'t mut FunctionBuilderContext) -> Self { Self{
			builder: frontend::FunctionBuilder::new(function, function_builder_context),
			constants_u32: default(),
			constants_f32: default(),
			//constants_f64: default(),
	} }
	fn u32(&mut self, value: u32) -> Value {
		match self.constants_u32.entry(value) {
			std::collections::hash_map::Entry::Occupied(value) => *value.get(),
			std::collections::hash_map::Entry::Vacant(entry) => *entry.insert(self.builder.ins().iconst(I32, value as i64))
		}
	}
	fn f32(&mut self, value: f32) -> Value {
		match self.constants_f32.entry(value.to_bits()) {
			std::collections::hash_map::Entry::Occupied(value) => *value.get(),
			std::collections::hash_map::Entry::Vacant(entry) => *entry.insert(self.builder.ins().f32const(value))
		}
	}
}

fn load(base: Value, index: usize, f: &mut FunctionBuilder) -> Value { f.ins().load(F32, MemFlags::trusted(), base, (index*std::mem::size_of::<f32>()) as i32) }
/*fn store(&mut self, value: Value, base: Value, index: usize) { self.ins().store(MemFlags::trusted(), value, base, (index*std::mem::size_of::<f64>()) as i32); }
fn neg(&mut self, x: Value) -> Value { self.ins().fneg(x) }
fn min(&mut self, a: Value, b: Value) -> Value { self.ins().fmin(a, b) }*/
fn max(a: Value, b: Value, f: &mut FunctionBuilder) -> Value { f.ins().fmax(a, b) }
fn add(a: Value, b: Value, f: &mut FunctionBuilder) -> Value { f.ins().fadd(a, b) }
fn sub(a: Value, b: Value, f: &mut FunctionBuilder) -> Value { f.ins().fsub(a, b) }
fn mul(a: Value, b: Value, f: &mut FunctionBuilder) -> Value { f.ins().fmul(a, b) }
fn div(a: Value, b: Value, f: &mut FunctionBuilder) -> Value { f.ins().fdiv(a, b) }
fn fma(a: Value, b: Value, c: Value, f: &mut FunctionBuilder) -> Value { f.ins().fma(a,b,c) } //let mul = self.mul(a, b); self.add(mul, c) }

macro_rules! f {
	[$f:ident $function:ident($arg0:expr)] => {{let arg0 = $arg0; $f.ins().$function(arg0)}};
	[$f:ident $function:ident($arg0:expr, $arg1:expr)] => {{let arg0 = $arg0; let arg1 = $arg1; $f.ins().$function(arg0, arg1)}};
}

fn exp2(x: Value, f: &mut FunctionBuilder) -> Value {
	//let x = f.ins().fdemote(F32, x);
	let x = max(x, f.f32(-126.99999), f);
	let ipart = f![f fcvt_to_sint(I32, sub(x, f.f32(1./2.), f))];
	let fpart = sub(x, f.ins().fcvt_from_sint(F32, ipart), f);
	let expipart = f![f bitcast(F32, f![f ishl_imm(f![f iadd(ipart, f.u32(127))], 23)])];
	let exp = [0.99999994f32, 0.69315308, 0.24015361, 0.055826318, 0.0089893397, 0.0018775767].map(|c| f.f32(c));
	let expfpart = fma(fma(fma(fma(fma(exp[5], fpart, exp[4], f), fpart, exp[3], f), fpart, exp[2], f), fpart, exp[1], f), fpart, exp[0], f);
	mul(expipart, expfpart, f)
	//f![f fpromote(F64, f.mul(expipart, expfpart))]
}
fn log2(x: Value, f: &mut FunctionBuilder) -> Value {
	//let x = f.ins().fdemote(F32, x);
	let i = f.ins().bitcast(I32, x);
	let e = f![f fcvt_from_sint(F32, f![f isub(f![f ushr_imm(f![f band(i, f.u32(0x7F800000))], 23)], f.u32(127))])];
	let m = f![f bitcast(F32, f![f bor(f![f band(i, f.u32(0x007FFFFF))], f.u32(1f32.to_bits()))])];
	let log = [3.1157899f32, -3.3241990, 2.5988452, -1.2315303,  0.31821337, -0.034436006].map(|c| f.f32(c));
	let p = fma(fma(fma(fma(fma(log[5], m, log[4], f), m, log[3], f), m, log[2], f), m, log[1], f), m, log[0], f);
	let p = mul(p, sub(m, f.f32(1.), f), f); //?
	//f![f fpromote(F64, f.add(p, e))]
	add(p, e, f)
}

#[derive(derive_more::Deref,derive_more::DerefMut)] struct AstFunctionBuilder<'t> {
	#[deref]#[deref_mut] builder: FunctionBuilder<'t>,
	definitions: linear_map::LinearMap<ast::Value, Value>,
}

use ast::*;
impl AstFunctionBuilder {
fn expr(&mut self, e: &Expression) -> Value {
	use Expression::*;
	match e {
		&Literal(v) => self.f32(v as f32),
		Use(v) => *self.definitions.get(v).unwrap_or_else(|| panic!("{v:?}")),
		Mul(a, b) => mul(self.expr(a), self.expr(b), self),
		Div(a, b) => div(self.expr(a), self.expr(b), self),
		Call { function, arguments } => match *function {
			"exp2" => exp2(self.expr(&arguments[0]), self),
			"log2" => log2(self.expr(&arguments[0]), self),
			function => panic!("{}", function)
		},
		_ => panic!("{e:?}")
	}
}
fn push(&mut self, statement: &Statement) {
	use Statement::*;
	match statement {
		Define { id, value } => { let value = expr(f, value); self.definitions.insert(ast::Value(id.0), value); },
		Branch { condition, consequent, alternative } => {
			let condition = self.expr(condition);
			let consequent_block = self.create_block();
			self.ins().brnz(condition, consequent_block, &[]);
			for s in &**consequent { f.push(s) }
			self.instructions.push("} else {".into());
			for s in &**alternative { f.push(s) }
			self.instructions.push("}".into());
		}
		_ => panic!("{statement:?}")
	}
}

pub fn compile(ast: &ast::Function) -> Function {
	let mut function = Function::new();
	let mut function_builder_context = FunctionBuilderContext::new();
	let mut f = FunctionBuilder::new(&mut function, &mut function_builder_context);
	let entry_block = f.create_block();
	f.func.signature.params = vec![AbiParam::new(I64); 2];
	f.append_block_params_for_function_params(entry_block);
	f.switch_to_block(entry_block);
	f.seal_block(entry_block);
	let_!{ [input, _output] = f.block_params(entry_block) => {
	let input = *input;
	let definitions = (0..ast.input).map(|i| (Value(i), load(input, i, &mut f))).collect();
	let ref mut f = AstFunctionBuilder{builder: f, definitions};
	for statement in &*ast.statements { f.push(statement); }
	f.ins().return_(&[]);
	function
}}}

#[cfg(feature="jit")] pub fn assemble(function: Function) -> impl Fn(&[f32], &mut [f32]) {
	let mut module = cranelift_jit::JITModule::new({
		let flag_builder = cranelift_codegen::settings::builder();
		use cranelift_codegen::settings::Configurable;
		let mut flag_builder = flag_builder;
		flag_builder.enable("is_pic").unwrap();
		flag_builder.set("enable_probestack", "false").unwrap();
		cranelift_jit::JITBuilder::with_isa(cranelift_native::builder().unwrap().finish(cranelift_codegen::settings::Flags::new(flag_builder)), cranelift_module::default_libcall_names())
	});
	use cranelift_module::Module;
	let mut context = module.make_context();
	context.func = function;
  let id = module.declare_function(&"", cranelift_module::Linkage::Export, &context.func.signature).unwrap();
  module.define_function(id, &mut context, &mut cranelift_codegen::binemit::NullTrapSink{}, &mut cranelift_codegen::binemit::NullStackMapSink{}).unwrap();
	module.finalize_definitions();
	let function = module.get_finalized_function(id);
	let function = unsafe{std::mem::transmute::<_,extern fn(*const f32, *mut f32)>(function)};
	move |input, output| function(input.as_ptr(), output.as_mut_ptr())
}
