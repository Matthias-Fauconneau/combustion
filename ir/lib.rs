//use std::default::default;
use ast::*;
pub use cranelift::codegen::ir::{function::Function};//, types::{I32, I64, F32, F64}};
/*use cranelift::{
	frontend::{self as frontend, FunctionBuilderContext},
	codegen::{ir::{types::Type, InstBuilder, MemFlags, entities::{Value, SigRef}, AbiParam, ExternalName, ExtFuncData}},
};

fn import(function: &mut Function, index: u32, signature: SigRef) {
	assert!(function.import_function(ExtFuncData{name: ExternalName::User{namespace: 0, index}, signature, colocated: default()}).as_u32() == index)
}

#[derive(derive_more::Deref,derive_more::DerefMut)] struct FunctionBuilder<'t> {
	#[deref]#[deref_mut] builder: frontend::FunctionBuilder<'t>,
	constants_u32: std::collections::HashMap<u32, Value>,
	constants_f32: std::collections::HashMap<u32, Value>,
	constants_f64: std::collections::HashMap<u64, Value>,
}

impl FunctionBuilder<'t> {
	fn new(function: &'t mut Function, function_builder_context: &'t mut FunctionBuilderContext) -> Self { Self{
			builder: frontend::FunctionBuilder::new(function, function_builder_context),
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
  /*fn load(&mut self, base: Value, index: usize) -> Value { self.ins().load(F64, MemFlags::trusted(), base, (index*std::mem::size_of::<f64>()) as i32) }
	fn store(&mut self, value: Value, base: Value, index: usize) { self.ins().store(MemFlags::trusted(), value, base, (index*std::mem::size_of::<f64>()) as i32); }
	fn neg(&mut self, x: Value) -> Value { self.ins().fneg(x) }
	fn min(&mut self, x: Value, y: Value) -> Value { self.ins().fmin(x, y) }
	fn max(&mut self, x: Value, y: Value) -> Value { self.ins().fmax(x, y) }
	fn add(&mut self, x: Value, y: Value) -> Value { self.ins().fadd(x, y) }
	fn sub(&mut self, x: Value, y: Value) -> Value { self.ins().fsub(x, y) }
	fn mul(&mut self, x: Value, y: Value) -> Value { self.ins().fmul(x, y) }
	fn div(&mut self, x: Value, y: Value) -> Value { self.ins().fdiv(x, y) }
	fn fma(&mut self, x: Value, y: Value, z: Value) -> Value { let mul = self.mul(x, y); self.add(mul, z) }*/
	fn exp2(&mut self, x: Value) -> Value {
		let f = self;
		let x = f.ins().fdemote(F32, x);
		let x = max(x, f.f32(-126.99999), f);
		let ipart = f![f fcvt_to_sint(I32, sub(x, f.f32(1./2.), f))];
		let fpart = sub(x, f.ins().fcvt_from_sint(F32, ipart), f);
		let expipart = f![f bitcast(F32, f![f ishl_imm(f![f iadd(ipart, f.u32(127))], 23)])];
		let exp = [0.99999994f32, 0.69315308, 0.24015361, 0.055826318, 0.0089893397, 0.0018775767].map(|c| f.f32(c));
		let expfpart = fma(fma(fma(fma(f.fma(exp[5], fpart, exp[4]), fpart, exp[3], f), fpart, exp[2], f), fpart, exp[1], f), fpart, exp[0], f);
		f![f fpromote(F64, f.mul(expipart, expfpart))]
	}
	fn log2(&mut self, x: Value) -> Value {
		let f = self;
		let x = f.ins().fdemote(F32, x);
		let i = f.ins().bitcast(I32, x);
		let e = f![f fcvt_from_sint(F32, f![f isub(f![f ushr_imm(f![f band(i, f.u32(0x7F800000))], 23)], f.u32(127))])];
		let m = f![f bitcast(F32, f![f bor(f![f band(i, f.u32(0x007FFFFF))], f.u32(1f32.to_bits()))])];
		let log = [3.1157899f32, -3.3241990, 2.5988452, -1.2315303,  0.31821337, -0.034436006].map(|c| f.f32(c));
		let p = fma(fma(fma(fma(f.fma(log[5], m, log[4]), m, log[3], f), m, log[2], f), m, log[1], f), m, log[0], f);
		let p = mul(p, sub(m, f.f32(1.), f), f); //?
		f![f fpromote(F64, f.add(p, e))]
	}
}*/

pub fn compile<const V: usize, const A: usize>(s: &Subroutine<V, A>) -> Function {
	let mut function = Function::new();
	let mut function_builder_context = FunctionBuilderContext::new();
	let ref mut f = FunctionBuilder::new(&mut function, &mut function_builder_context);
	let entry_block = f.create_block();
	const parameters: [Type; {V+A}] = [[F64; V], [I64; A]];
	f.func.signature.params = map(&parameters, |t| AbiParam::new(*t)).to_vec();
	f.append_block_params_for_function_params(entry_block);
	f.switch_to_block(entry_block);
	f.seal_block(entry_block);
	let parameters = f.block_params(entry_block);
	for statement in statements {
		use Statement::*;
		match statement {
			/*ConditionalStatement { condition, consequent, alternative } => {
				if self.eval(condition) != 0. { self.run(consequent); } else { self.run(alternative); }
			},
			Definition { id, value } => {
				let value = self.eval(value);
				assert!(self.definitions.insert(Value(id.0), value).is_none());
			},
			Output { index, value } => self.output[*index] = self.eval(value),*/
			_ => unimplemented!()
		}
	}
	f.ins().return_(&[]);
	function
}

#[cfg(feature="jit")] pub fn assemble<const V: usize, const A: usize>(function: Function) -> impl Fn([f64; V], [&[f64]; A], &mut [f64]) {
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
	let function = unsafe{std::mem::transmute::<_,extern fn(*const f64/*[V]*/, *const *const f64/*[A]*/, *mut f64)>(function)};
	move |value, array, output| function(value.as_ptr(), array.as_ptr(), output.as_mut_ptr());
}
