#![feature(default_free_fn,format_args_capture,in_band_lifetimes,array_map)]
use {std::default::default, iter::map};
pub use cranelift::codegen::ir::{function::Function, types::{Type, I32, I64, F32, F64}, condcodes::FloatCC};
use cranelift::{
	frontend::{FunctionBuilder, FunctionBuilderContext},
	codegen::{ir::{InstBuilder, AbiParam, entities::Value, MemFlags}, }, //types::Type,
};

#[derive(derive_more::Deref,derive_more::DerefMut)] struct Builder<'t> {
	#[deref]#[deref_mut] builder: FunctionBuilder<'t>,
	constants_u32: std::collections::HashMap<u32, Value>,
	constants_f32: std::collections::HashMap<u32, Value>,
	constants_f64: std::collections::HashMap<u64, Value>,
}

impl Builder<'t> {
	fn new(function: &'t mut Function, function_builder_context: &'t mut FunctionBuilderContext) -> Self { Self{
			builder: FunctionBuilder::new(function, function_builder_context),
			constants_u32: default(),
			constants_f32: default(),
			constants_f64: default(),
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
	fn f64(&mut self, value: f64) -> Value {
		match self.constants_f64.entry(value.to_bits()) {
			std::collections::hash_map::Entry::Occupied(value) => *value.get(),
			std::collections::hash_map::Entry::Vacant(entry) => *entry.insert(self.builder.ins().f64const(value))
		}
	}
}

fn load(base: Value, index: usize, f: &mut Builder) -> Value { f.ins().load(F64, MemFlags::trusted(), base, (index*std::mem::size_of::<f64>()) as i32) }
//fn cast(to: Type, x: Value, f: &mut Builder) -> Value { f.ins().bitcast(to, x) }
fn and(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().band(a, b) }
fn or(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().bor(a, b) }
fn ishl_imm(x: Value, imm: u8, f: &mut Builder) -> Value { f.ins().ishl_imm(x, imm as i64) }
fn ushr_imm(x: Value, imm: u8, f: &mut Builder) -> Value { f.ins().ushr_imm(x, imm as i64) }
fn iadd(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().iadd(a, b) }
fn isub(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().isub(a, b) }
fn neg(x: Value, f: &mut Builder) -> Value { f.ins().fneg(x) }
//fn min(&mut self, a: Value, b: Value) -> Value { self.ins().fmin(a, b) }
fn max(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fmax(a, b) }
fn add(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fadd(a, b) }
fn sub(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fsub(a, b) }
fn mul(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fmul(a, b) }
//fn fma(a: Value, b: Value, c: Value, f: &mut Builder) -> Value { f.ins().fma(a,b,c) }
fn fma(a: Value, b: Value, c: Value, f: &mut Builder) -> Value { add(mul(a, b, f), c, f) }
fn div(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fdiv(a, b) }
fn sqrt(x: Value, f: &mut Builder) -> Value { f.ins().sqrt(x) }
fn fpromote(x: Value, f: &mut Builder) -> Value { f.ins().fpromote(F64, x) }
fn fdemote(x: Value, f: &mut Builder) -> Value { f.ins().fdemote(F32, x) }
fn fcvt_to_sint(x: Value, f: &mut Builder) -> Value { f.ins().fcvt_to_sint(I32, x) }
fn fcvt_from_sint(x: Value, f: &mut Builder) -> Value { f.ins().fcvt_from_sint(F32, x) }
fn store(value: Value, base: Value, index: usize, f: &mut Builder) { f.ins().store(MemFlags::trusted(), value, base, (index*std::mem::size_of::<f64>()) as i32); }

#[derive(derive_more::Deref,derive_more::DerefMut)] struct AstBuilder<'t> {
	#[deref]#[deref_mut] builder: Builder<'t>,
	types: linear_map::LinearMap<ast::Value, ast::Type>,
	values: linear_map::LinearMap<ast::Value, Value>,
}

use ast::*;
impl AstBuilder<'_> {
fn pass(&mut self, e: &Expression) -> ast::Type { // check_types_and_load_constants
	use Expression::*;
	match e {
		&F32(v) => { self.f32(v as f32); ast::Type::F32 },
		&F64(v) => { self.f64(v); ast::Type::F64 },
		&I32(v) => { self.u32(v); ast::Type::I32 },
		Value(v) => *self.types.get(v).unwrap_or_else(|| panic!("{:?} {v:?}", self.types)),
		Cast(to, x) => { self.pass(x); *to },
		Neg(x)|IShLImm(x,_)|UShRImm(x,_)|Sqrt(x) => self.pass(x),
		FPromote(x) => { self.pass(x); ast::Type::F64 },
		FDemote(x) => { self.pass(x); ast::Type::F32 },
		FCvtToSInt(x)  => { self.pass(x); ast::Type::I32 }
		FCvtFromSInt(x) => { self.pass(x); ast::Type::F32 },
		And(a,b)|Or(a,b)|IAdd(a,b)|ISub(a,b)|Max(a,b)|Add(a,b)|Sub(a,b)|Mul(a,b)|Div(a,b)|LessOrEqual(a,b) => { let [a,b] = [a,b].map(|x| self.pass(x)); assert!(a==b,"{e:?}"); a },
		MulAdd(a,b,c) => { let [a,b,c] = [a,b,c].map(|x| self.pass(x)); assert!(a==b && b==c); a },
		Block { statements, result } => { for s in &**statements { self.check_types_and_load_constants(s) } self.pass(result) }
	}
}
fn check_types_and_load_constants(&mut self, statement: &Statement) {
	use Statement::*;
	match statement {
		Value { id, value } => { let value = self.load(value); assert!(!self.types.contains_key(id)); self.types.insert(id.clone(), value); },
		Branch { condition, true_exprs, false_exprs, results } => {
			self.load(condition);
			// Always load constants of both branches in root block so that they are always available for all later uses
			let true_types = map(&**true_exprs, |e| self.pass(e));
			let false_types = map(&**false_types, |e| self.pass(e));
			for id in &**results { assert!(!self.values.contains_key(id)); }
			self.types.extend(results.iter().zip(true_types.iter().zip(&*false_types)).map(|(id, (&a,&b))| { assert!(a==b); (id.clone(), a) }));
		}
	}
}
fn expr(&mut self, e: &Expression) -> Value {
	use Expression::*;
	match e {
		Value(v) => *self.values.get(v).unwrap_or_else(|| panic!("{:?} {v:?}", self.values)),
		&F32(v) => self.constants_f32[&(v as f32).to_bits()],
		&F64(v) => self.constants_f64[&v.to_bits()],
		I32(v) => self.constants_u32[v],
		Cast(to, x) => unimplemented!("{to:?} {x:?}"),//cast(match to {ast::Type::F32=>self::F32,ast::Type::F64=>self::F64,ast::Type::I32=>I32}, self.expr(x), self),
		And(a, b) => and(self.expr(a), self.expr(b), self),
		Or(a, b)  => or(self.expr(a), self.expr(b), self),
		IShLImm(x, imm) => ishl_imm(self.expr(x), *imm, self),
		UShRImm(x, imm) => ushr_imm(self.expr(x), *imm, self),
		IAdd(a, b) => iadd(self.expr(a), self.expr(b), self),
		ISub(a, b) => isub(self.expr(a), self.expr(b), self),
		Neg(x) => neg(self.expr(x), self),
		Max(a, b) => max(self.expr(a), self.expr(b), self),
		LessOrEqual(a, b) => { let (a, b) = (self.expr(a), self.expr(b)); self.ins().fcmp(FloatCC::LessThanOrEqual, a, b) },
		Add(a, b) => add(self.expr(a), self.expr(b), self),
		Sub(a, b) => sub(self.expr(a), self.expr(b), self),
		Mul(a, b) => mul(self.expr(a), self.expr(b), self),
		MulAdd(a, b, c) => fma(self.expr(a), self.expr(b), self.expr(c), self),
		Div(a, b) => div(self.expr(a), self.expr(b), self),
		Sqrt(x) => sqrt(self.expr(x), self),
		FPromote(x) => fpromote(self.expr(x), self),
		FDemote(x) => fdemote(self.expr(x), self),
		FCvtToSInt(x) => fcvt_to_sint(self.expr(x), self),
		FCvtFromSInt(x) => fcvt_from_sint(self.expr(x), self),
		Block { statements, result } => {
			for s in &**statements { self.push(s) }
			self.expr(result)
		}
	}
}
fn push(&mut self, statement: &Statement) {
	use Statement::*;
	match statement {
		Value { id, value } => { let value = self.expr(value); assert!(!self.values.contains_key(id)); self.values.insert(id.clone(), value); },
		Branch { condition, true_exprs, false_exprs, results } => {
			let condition = self.expr(condition);
			/*// Always load constants of both branches in root block so that they are always available for all later uses (already done for non-root branches (in check_types_and_load_constants))
			for e in &**true_exprs { self.load(e) }
			for e in &**false_exprs { self.load(e) }*/
			let false_block = self.create_block();
			self.ins().brz(condition, false_block, &[]);
			self.seal_block(false_block);
			let true_block = self.create_block();
			self.ins().jump(true_block, &[]);
			self.seal_block(true_block);
			let join = self.create_block();
			for _ in 0..results.len() { self.append_block_param(join, F64); }

			let scope = self.values.clone();
			self.switch_to_block(true_block);
			let true_values = map(&**true_exprs, |e| self.expr(e));
			self.values = scope;
			assert!(results.len() == true_values.len());
			self.ins().jump(join, &true_block);

			let scope = self.values.clone();
			self.switch_to_block(false_block);
			let false_values = map(&**false_exprs, |e| self.expr(e));
			self.values = scope;
			assert!(results.len() == false_values.len());
			self.ins().jump(join, &false_block);

			self.seal_block(join);
			self.switch_to_block(join);
			let params = self.builder.block_params(join);
			assert!(results.len() == params.len());
			for id in &**results { assert!(!self.values.contains_key(id)); }
			self.values.extend(results.iter().zip(params).map(|(id, &value)| (id.clone(), value)));
		}
	}
}
}

pub fn compile(ast: &ast::Function) -> Function {
	let mut function = Function::new();
	let mut function_builder_context = FunctionBuilderContext::new();
	let mut f = Builder::new(&mut function, &mut function_builder_context);
	let entry_block = f.create_block();
	f.func.signature.params = vec![AbiParam::new(I64); 2];
	f.append_block_params_for_function_params(entry_block);
	f.switch_to_block(entry_block);
	f.seal_block(entry_block);
	let_!{ &[input, output] = f.block_params(entry_block) => {
	let types = (0..ast.input).map(|i| (Value(i), ast::Type::F64)).collect();
	let values = (0..ast.input).map(|i| (Value(i), load(input, i, &mut f))).collect();
	let ref mut f = AstBuilder{builder: f, types, values};
	for statement in &*ast.statements { f.check_types_and_load_constants(statement); }
	for statement in &*ast.statements { f.push(statement); }
	for (i, e) in ast.output.iter().enumerate() { store(f.expr(e), output, i, f); }
	f.ins().return_(&[]);
	function
}}}

#[cfg(feature="jit")] pub fn assemble(function: Function) -> impl Fn(&[f64], &mut [f64]) {
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
	let function = unsafe{std::mem::transmute::<_,extern fn(*const f64, *mut f64)>(function)};
	move |input, output| function(input.as_ptr(), output.as_mut_ptr())
}
