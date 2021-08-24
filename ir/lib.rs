#![feature(default_free_fn,format_args_capture)]
use {std::default::default, iter::map};
pub use cranelift::codegen::ir::{function::Function, types::{Type, I32, I64, F32, F64}, condcodes::FloatCC};
use cranelift::{
	frontend::{FunctionBuilder, FunctionBuilderContext},
	codegen::{ir::{InstBuilder, AbiParam, entities::Value, MemFlags}, }, //types::Type,
};

use ordered_float::NotNan;
type R32 = NotNan<f32>;
type R64 = NotNan<f64>;

#[derive(derive_more::Deref,derive_more::DerefMut)] struct Builder<'t> {
	#[deref]#[deref_mut] builder: FunctionBuilder<'t>,
	//constants_u32: std::collections::HashMap<u32, Value>,
	constants_f32: std::collections::HashMap<R32, Value>,
	constants_f64: std::collections::HashMap<R64, Value>,
}

impl<'t> Builder<'t> {
	fn new(function: &'t mut Function, function_builder_context: &'t mut FunctionBuilderContext) -> Self { Self{
			builder: FunctionBuilder::new(function, function_builder_context),
			//constants_u32: default(),
			constants_f32: default(),
			constants_f64: default(),
	} }
	//fn u32(&mut self, value: u32) -> Value { *self.constants_u32.entry(value).or_insert_with(|| self.builder.ins().iconst(I32, value as i64)) }
	fn f32(&mut self, value: f32) -> Value { *self.constants_f32.entry(R32::new(value).unwrap()).or_insert_with(|| self.builder.ins().f32const(value)) }
	fn f64(&mut self, value: f64) -> Value { *self.constants_f64.entry(R64::new(value).unwrap()).or_insert_with(|| self.builder.ins().f64const(value)) }
}

fn load(r#type: Type, base: Value, index: usize, f: &mut Builder) -> Value { f.ins().load(r#type, MemFlags::trusted(), base, (index as u32*r#type.bytes()) as i32) }
//fn cast(to: Type, x: Value, f: &mut Builder) -> Value { f.ins().bitcast(to, x) }
/*fn and(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().band(a, b) }
fn or(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().bor(a, b) }
fn ishl_imm(x: Value, imm: u8, f: &mut Builder) -> Value { f.ins().ishl_imm(x, imm as i64) }
fn ushr_imm(x: Value, imm: u8, f: &mut Builder) -> Value { f.ins().ushr_imm(x, imm as i64) }
fn iadd(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().iadd(a, b) }
fn isub(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().isub(a, b) }*/
fn neg(x: Value, f: &mut Builder) -> Value { f.ins().fneg(x) }
//fn min(&mut self, a: Value, b: Value) -> Value { self.ins().fmin(a, b) }
fn max(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fmax(a, b) }
fn add(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fadd(a, b) }
fn sub(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fsub(a, b) }
fn mul(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fmul(a, b) }
//fn fma(a: Value, b: Value, c: Value, f: &mut Builder) -> Value { f.ins().fma(a,b,c) }
//fn fma(a: Value, b: Value, c: Value, f: &mut Builder) -> Value { add(mul(a, b, f), c, f) }
fn div(a: Value, b: Value, f: &mut Builder) -> Value { f.ins().fdiv(a, b) }
fn sqrt(x: Value, f: &mut Builder) -> Value { f.ins().sqrt(x) }
/*fn fpromote(x: Value, f: &mut Builder) -> Value { f.ins().fpromote(F64, x) }
fn fdemote(x: Value, f: &mut Builder) -> Value { f.ins().fdemote(F32, x) }
fn fcvt_to_sint(x: Value, f: &mut Builder) -> Value { f.ins().fcvt_to_sint(I32, x) }
fn fcvt_from_sint(x: Value, f: &mut Builder) -> Value { f.ins().fcvt_from_sint(F32, x) }*/
fn store(value: Value, base: Value, index: usize, f: &mut Builder) { f.ins().store(MemFlags::trusted(), value, base, (index*std::mem::size_of::<f32>()) as i32); }

fn exp_approx_constants(f: &mut Builder) -> Type { //e-12 (19*,1/) (-9->-7)
	f.f32(1./2048.);f.f32(1.);f.f32(3./28.);f.f32(1./1680.);f.f32(1./2.);f.f32(1./84.);
	F32
}
fn exp_approx(x: Value, f: &mut Builder) -> Value { //e-12 (19*,1/) (-9->-7)
	let x = mul(f.f32(1./2048.), x, f);
	let x2 = mul(x, x, f);
	let x3 = mul(x, x2, f);
	let a = add(f.f32(1.), add(mul(f.f32(3./28.), x2, f), mul(f.f32(1./1680.), mul(x, x3, f),f),f),f);
	let b = add(mul(f.f32(1./2.), x, f), mul(f.f32(1./84.), x3, f), f);
	let sq = |x,f:&mut Builder| { mul(x,x,f) };
	sq(sq(sq(sq(sq(sq(sq(sq(sq(sq(sq(div(add(a,b,f),sub(a,b,f),f),f),f),f),f),f),f),f),f),f),f),f)
}

fn ln_approx_constants(x0: f64, f: &mut Builder) -> Type {
	f.f32(1.);f.f32(f64::ln(x0) as _);f.f32(16.*2.);f.f32(1./3.);f.f32(1./5.);f.f32(1./7.);f.f32(1./9.);
	F32
}
fn ln_approx(x0: f64, x: Value, f: &mut Builder) -> Value { // -5
	let x = mul(f.f32((1./x0) as _), x, f);
	let x = sqrt(sqrt(sqrt(sqrt(x,f),f),f),f);
	let x = div(sub(x,f.f32(1.),f),add(x,f.f32(1.),f),f);
	let x2 = mul(x,x,f);
	let x4 = mul(x2,x2,f);
	let x6 = mul(x4,x2,f);
	add(f.f32(f64::ln(x0) as _), mul(mul(f.f32(16.*2.),x,f), add(f.f32(1.), add(add(add(mul(f.f32(1./3.), x2, f), mul(f.f32(1./5.),x4,f), f), mul(f.f32(1./7.),x6,f),f), mul(f.f32(1./9.),mul(x6,x2,f),f),f),f),f),f)
}

#[derive(derive_more::Deref,derive_more::DerefMut)] struct AstBuilder<'t,'m> {
	#[deref]#[deref_mut] builder: &'m mut Builder<'t>,
	types: linear_map::LinearMap<ast::Value, ast::Type>,
	values: linear_map::LinearMap<ast::Value, Value>,
}

use ast::*;
impl AstBuilder<'_,'_> {
/*fn inline_pass(&mut self, x: &Expression, function: impl Fn(Expression, &mut Block) -> Expression) -> ast::Type {
	self.pass(x);
	let ref mut names = vec!["x".to_string()];
	let ref mut block = Block{names, statements: vec![]};
	let result = function(Expr::Value(ast::Value(0)).into(), block),
	for statement in &*block.statements { self.check_types_and_load_constants(statement); }
	f.pass(&result)
}
fn inline_expr(&mut self, x: &Expression, function: impl Fn(Expression, &mut Block) -> Expression) -> Value {
	let types = vec![(ast::Value(0), float::TYPE)].into_iter().collect();
	let values = vec![(ast::Value(0), self.expr(x))].into_iter().collect();
	let f = AstBuilder{builder: &mut self.builder, types, values}; // New values scope/frame to let function define new intermediates without conflicting with parent caller
	let ref mut names = vec!["x".to_string()];
	let ref mut block = Block{names, statements: vec![]};
	let result = function(Expr::Value(ast::Value(0)).into(), block);
	for statement in &*block.statements { f.push(statement); }
	f.expr(&result)
}*/
fn pass(&mut self, e: &Expression) -> ast::Type { // check_types_and_load_constants
	use Expr::*;
	match e { Expression::Expr(e) => match e {
		//&I32(v) => { self.u32(v); ast::Type::I32 },
		&F32(v) => { self.f32(*v); ast::Type::F32 },
		&F64(v) => { self.f64(*v); ast::Type::F64 },
		Value(v) => *self.types.get(v).unwrap_or_else(|| panic!("{:?} {v:?}", self.types)),
		//Cast(to, x) => { self.pass(x); *to },
		Neg(x)/*|IShLImm(x,_)|UShRImm(x,_)*/|Sqrt(x)|Sq(x) => self.pass(x),
		Exp(x) => {self.pass(x); exp_approx_constants(self); ast::Type::F32}, //self.inline_pass(x, exp_approx),
		Ln{x0,x} => {self.pass(x); ln_approx_constants(**x0, self); ast::Type::F32},  //self.inline_pass(x, |x,f| ln_approx(**x0, x, f)),
		/*FPromote(x) => { self.pass(x); ast::Type::F64 },
		FDemote(x) => { self.pass(x); ast::Type::F32 },
		FCvtToSInt(x)  => { self.pass(x); ast::Type::I32 }
		FCvtFromSInt(x) => { self.pass(x); ast::Type::F32 },*/
		/*And(a,b)|Or(a,b)|IAdd(a,b)|ISub(a,b)|*/Max(a,b)|Add(a,b)|Sub(a,b)|Mul(a,b)|Div(a,b)|LessOrEqual(a,b) => { let [a,b] = [a,b].map(|x| self.pass(x)); assert!(a==b,"{e:?}"); a },
		//MulAdd(a,b,c) => { let [a,b,c] = [a,b,c].map(|x| self.pass(x)); assert!(a==b && b==c); a },
		},
		Expression::Block { statements, result } => { for s in &**statements { self.check_types_and_load_constants(s) } self.pass(result) }
	}
}
fn check_types_and_load_constants(&mut self, statement: &Statement) {
	use Statement::*;
	match statement {
		Value { id, value } => { let value = self.pass(value); assert!(!self.types.contains_key(id)); self.types.insert(id.clone(), value); },
		Select { condition, true_exprs, false_exprs, results } => {
			self.pass(condition);
			// Always load constants of both branches in root block so that they are always available for all later uses
			let true_types = map(&**true_exprs, |e| self.pass(e));
			let false_types = map(&**false_exprs, |e| self.pass(e));
			for id in &**results { assert!(!self.values.contains_key(id)); }
			self.types.extend(results.iter().zip(true_types.iter().zip(&*false_types)).map(|(id, (&a,&b))| { assert!(a==b); (id.clone(), a) }));
		}
		//Display(_) => unimplemented!(),
	}
}
fn expr(&mut self, e: &Expression) -> Value {
	use Expr::*;
	match e { Expression::Expr(e) => match e {
		Value(v) => *self.values.get(v).unwrap_or_else(|| panic!("{:?} {v:?}", self.values)),
		F32(v) => self.constants_f32[v],
		F64(v) => self.constants_f64[&R64::new(**v).unwrap()],
		/*I32(v) => self.constants_u32[v],
		Cast(to, x) => unimplemented!("{to:?} {x:?}"),//cast(match to {ast::Type::F32=>self::F32,ast::Type::F64=>self::F64,ast::Type::I32=>I32}, self.expr(x), self),
		And(a, b) => and(self.expr(a), self.expr(b), self),
		Or(a, b)  => or(self.expr(a), self.expr(b), self),
		IShLImm(x, imm) => ishl_imm(self.expr(x), *imm, self),
		UShRImm(x, imm) => ushr_imm(self.expr(x), *imm, self),
		IAdd(a, b) => iadd(self.expr(a), self.expr(b), self),
		ISub(a, b) => isub(self.expr(a), self.expr(b), self),*/
		Neg(x) => neg(self.expr(x), self),
		Max(a, b) => max(self.expr(a), self.expr(b), self),
		LessOrEqual(a, b) => { let (a, b) = (self.expr(a), self.expr(b)); self.ins().fcmp(FloatCC::LessThanOrEqual, a, b) },
		Add(a, b) => add(self.expr(a), self.expr(b), self),
		Sub(a, b) => sub(self.expr(a), self.expr(b), self),
		Mul(a, b) => mul(self.expr(a), self.expr(b), self),
		//MulAdd(a, b, c) => fma(self.expr(a), self.expr(b), self.expr(c), self),
		Div(a, b) => div(self.expr(a), self.expr(b), self),
		Sqrt(x) => sqrt(self.expr(x), self),
		/*FPromote(x) => fpromote(self.expr(x), self),
		FDemote(x) => fdemote(self.expr(x), self),
		FCvtToSInt(x) => fcvt_to_sint(self.expr(x), self),
		FCvtFromSInt(x) => fcvt_from_sint(self.expr(x), self),*/
		Exp(x) => exp_approx(self.expr(x), self), //self.inline_expr(x, exp_approx),
		Ln{x0,x} => ln_approx(**x0, self.expr(x), self), //self.inline_expr(x, |x,f| ln_approx(**x0, x, f)),
		Sq(x) => { let x = self.expr(x); mul(x, x, self)}, //self.inline_expr(x, |x,f| ln_approx(**x0, x, f)),
		}
		Expression::Block { statements, result } => {
			for s in &**statements { self.push(s) }
			self.expr(result)
		},
	}
}
fn push(&mut self, statement: &Statement) {
	use Statement::*;
	match statement {
		Value { id, value } => { let value = self.expr(value); assert!(!self.values.contains_key(id)); self.values.insert(id.clone(), value); },
		Select { condition, true_exprs, false_exprs, results } => {
			let condition = self.expr(condition);
			let false_block = self.create_block();
			self.ins().brz(condition, false_block, &[]);
			self.seal_block(false_block);
			let true_block = self.create_block();
			self.ins().jump(true_block, &[]);
			self.seal_block(true_block);
			let join = self.create_block();
			for _ in 0..results.len() { self.append_block_param(join, F32); }

			//let scope = self.values.clone();
			self.switch_to_block(true_block);
			let true_values = map(&**true_exprs, |e| self.expr(e));
			//self.values = scope;
			assert!(results.len() == true_values.len());
			self.ins().jump(join, &true_values);

			//let scope = self.values.clone();
			self.switch_to_block(false_block);
			let false_values = map(&**false_exprs, |e| self.expr(e));
			//self.values = scope;
			assert!(results.len() == false_values.len());
			self.ins().jump(join, &false_values);

			self.seal_block(join);
			self.switch_to_block(join);
			let params = self.builder.block_params(join);
			assert!(results.len() == params.len());
			for id in &**results { assert!(!self.values.contains_key(id)); }
			self.values.extend(results.iter().zip(params).map(|(id, &value)| (id.clone(), value)));
		}
		//Display(_) => unimplemented!()
	}
}
}

pub fn compile(ast: &ast::Function) -> Function {
	let mut function = Function::new();
	let mut function_builder_context = FunctionBuilderContext::new();
	let ref mut f = Builder::new(&mut function, &mut function_builder_context);
	let entry_block = f.create_block();
	f.func.signature.params = vec![AbiParam::new(I64); 2];
	f.append_block_params_for_function_params(entry_block);
	f.switch_to_block(entry_block);
	f.seal_block(entry_block);
	let_!{ &[input, output] = f.block_params(entry_block) => {
	let types = ast.input.iter().enumerate().map(|(i, &r#type)| (Value(i), r#type)).collect();
	let values = ast.input.iter().enumerate().map(|(i, r#type)| (Value(i), load({use ast::Type::*; match r#type {F32=>self::F32, F64=>self::F64}}, input, i, f))).collect();
	let ref mut f = AstBuilder{builder: f, types, values};
	for statement in &*ast.statements { f.check_types_and_load_constants(statement); }
	for statement in &*ast.statements { f.push(statement); }
	for (i, e) in ast.output.iter().enumerate() { f.pass(e); store(f.expr(e), output, i, f); }
	f.ins().return_(&[]);
	function
}}}

#[cfg(feature="jit")] pub fn assemble<T>(function: Function) -> impl Fn(&[T], &mut [T]) {
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
	let function = unsafe{std::mem::transmute::<_,extern fn(*const T, *mut T)>(function)};
	move |input, output| function(input.as_ptr(), output.as_mut_ptr())
}
