use super::*;

#[derive(PartialEq, Eq)] pub enum Property { Pressure, Volume }
#[derive(Clone, Copy)] pub struct Constant<const CONSTANT: Property>(pub f64);
#[derive(derive_more::Deref, Default)] pub struct StateVector<const CONSTANT: Property>(pub Box<[f64/*; T,P|V,[S-1]*/]>);
pub type Derivative<const CONSTANT: Property> = StateVector<CONSTANT>;

impl State {
	//pub fn constant<const CONSTANT: Property>(&Self{pressure, volume, ..}: &Self) -> f64 { // arbitrary_self_types
	pub fn constant<const CONSTANT: Property>(&self) -> Constant<CONSTANT> { let Self{pressure, volume, ..} = self;
		Constant(*{use Property::*; match CONSTANT {Pressure => pressure, Volume => volume}})
	}
}

use std::array::IntoIter;

impl<const CONSTANT: Property> From<&State> for StateVector<CONSTANT> {
	fn from(State{temperature, pressure, volume, amounts}: &State) -> Self {
		Self([*temperature, *{use Property::*; match CONSTANT {Pressure => volume, Volume => pressure}}].iter().chain(amounts[..amounts.len()-1].iter()).copied().collect())
	}
}

impl State {
	pub fn new<const CONSTANT: Property>(total_amount: f64, Constant(thermodynamic_state_constant): Constant<CONSTANT>, u: &StateVector<CONSTANT>) -> Self {
		let u = &u.0;
		assert!(!u.iter().any(|&n| n<0.));
		let amounts = &u[2..];
		let (pressure, volume) = {use Property::*; match CONSTANT {
			Pressure => (thermodynamic_state_constant, u[1]),
			Volume => (u[1], thermodynamic_state_constant)
		}};
		State{
			temperature: u[0],
			pressure: pressure,
			volume: volume,
			amounts: amounts.iter().chain(&[total_amount - iter::into::Sum::<f64>::sum(amounts)]).copied().collect()
		}
	}
}

impl std::fmt::Display for State {
	fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
		let Self{temperature, pressure, volume, amounts} = self;
		write!(fmt, "T: {}, P: {}, V: {}, n: {:?}", temperature/*/K*/, pressure*NA, volume, amounts)
	}
}

use {cranelift::prelude::{*, types::{/*I32,*/ F64}, codegen::{ir, binemit}}, cranelift_module::{Linkage, Module}};

macro_rules! f {
	[$f:ident $function:ident($arg0:expr)] => {{
		let arg0 = $arg0;
		$f.ins().$function(arg0)
	}};
	[$f:ident $function:ident($arg0:expr, $arg1:expr)] => {{
		let arg0 = $arg0;
		let arg1 = $arg1;
		$f.ins().$function(arg0, arg1)
	}};
	[$f:ident $function:ident($arg0:expr, $arg1:expr, $arg2:expr)] => {{
		let arg0 = $arg0;
		let arg1 = $arg1;
		let arg2 = $arg2;
		$f.ins().$function(arg0, arg1, arg2)
	}};
	[$f:ident $function:ident($arg0:expr, $arg1:expr, $arg2:expr, $arg3:expr)] => {{
		let arg0 = $arg0;
		let arg1 = $arg1;
		let arg2 = $arg2;
		let arg3 = $arg3;
		$f.ins().$function(arg0, arg1, arg2, arg3)
	}}
}

macro_rules! fma {
	[$f:ident ($arg0:expr, $arg1:expr, $arg2:expr)] => {{
		let arg0 = $arg0;
		let arg1 = $arg1;
		let mul = $f.ins().fmul(arg0, arg1);
		let arg2 = $arg2;
		$f.ins().fadd(mul, arg2)
	}};
}

/*trait Type { const r#type: types::Type; }
impl Type for i32 { const r#type : types::Type = types::I32; }
impl Type for f32 { const r#type : types::Type = types::F32; }
impl Type for f64 { const r#type : types::Type = types::F64; }
fn constant<T:/*std::fmt::Display+*/Type>(module: &mut cranelift_jit::JITModule, f: &mut FunctionBuilder<'_>, constant: T) -> Value {
where [(); std::mem::size_of::<T>()]: {
	let data = module.declare_data(&format!("{}: {}", constant, T::r#type), Linkage::Local, false, false).unwrap();
	module.define_data(
		data,
		&{let mut data_ctx = cranelift_module::DataContext::new(); data_ctx.define(Box::new(unsafe{std::mem::transmute_copy::<_,[u8; std::mem::size_of::<T>()]>(&constant)})); data_ctx}
	).unwrap();
	let data = module.declare_data_in_func(data, f.func);
	match T::r#type {
		r#type@types::F64 => f![f bitcast(r#type, f![f symbol_value(types::I64, data)])],  // FIXME: symbol_value(F64) ?
		r#type@types::F32 => f![f bitcast(r#type, f![f ireduce(types::I32, f![f symbol_value(types::I64, data)])])],  // FIXME: symbol_value(F32) ?
		r#type@types::I32 => f![f ireduce(r#type, f![f symbol_value(types::I64, data)])],  // FIXME: symbol_value(I32) ?
		r#type => f![f symbol_value(r#type, data)]
	} // Static data does not seem to work with JIT module
}*/
trait Load { fn load(&self, f: &mut FunctionBuilder) -> Value; }
/*impl Load for u32 { fn load(self, f: &mut FunctionBuilder<'_>) -> Value { f![f iconst(I32, self as i64)] } }
impl Load for i32 { fn load(self, f: &mut FunctionBuilder<'_>) -> Value { f![f iconst(I32, self as i64)] } }
impl Load for f32 { fn load(self, f: &mut FunctionBuilder<'_>) -> Value { f![f f32const(self)] } }*/
//impl Load for u32 { fn load(self, f: &mut FunctionBuilder<'_>) -> Value { f![f iconst(I32, self as i64)] } }
//impl Load for i32 { fn load(self, f: &mut FunctionBuilder<'_>) -> Value { f![f iconst(I32, self as i64)] } }
//impl Load for f32 { fn load(self, f: &mut FunctionBuilder<'_>) -> Value { f![f f32const(self)] } }
impl Load for f64 { fn load(&self, f: &mut FunctionBuilder<'_>) -> Value {
	//f![f vconst(F64, f.func.dfg.constants.insert(unsafe{std::mem::transmute_copy::<_,[u8; std::mem::size_of::<T>()]>(self).as_slice().into()}))]
	f![f f64const(*self)]
}}

struct Constants {
	//module: &'t mut cranelift_jit::JITModule,
	constants: std::collections::HashMap<u64, Value>,
	_1: Value,
	/*_1_f32: Value,
	_1_2: Value,
	_127: Value,
	exponent: Value,
	mantissa: Value,
	//exp: [Value; 3],
	//log: [Value; 3],
	exp: [Value; 6],
	log: [Value; 6],
	_m126_99999: Value,*/
}
impl Constants/*<'t>*/ {
	fn new(f: &mut FunctionBuilder/*<'t>*/) -> Self { Self{
		constants: Default::default(),
		_1: 1f64.load(f),
		/*_1_f32: 1f32.load(module, f),
		_1_2: (1./2f32).load(module, f),
		_127: 127.load(module, f),
		exponent: 0x7F800000u32.load(module, f),
		mantissa: 0x007FFFFFu32.load(module, f),
		//exp: [1.0017247, 0.65763628, 0.33718944].map(|c| f![f f32const(c)]),
		exp: [0.99999994f32, 0.69315308, 0.24015361, 0.055826318, 0.0089893397, 0.0018775767].map(|c| c.load(module, f)),
		//log: [2.28330284476918490682, -1.04913055217340124191, 0.204446009836232697516].map(|c| f![f f32const(c)]),
		log: [3.1157899f32, -3.3241990, 2.5988452, -1.2315303,  0.31821337, -0.034436006].map(|c| c.load(module, f)),
		_m126_99999: (-126.99999f32).load(module, f),*/
	}}
	#[track_caller] fn c(&mut self, f: &mut FunctionBuilder<'t>, constant: f64) -> Value {
		assert!(constant.is_finite(), "{}", constant);
		match self.constants.entry(constant.to_bits()) {
			std::collections::hash_map::Entry::Occupied(value) => *value.get(),
			std::collections::hash_map::Entry::Vacant(entry) => *entry.insert(constant.load(f))
		}
	}
}

pub extern fn _exp2(x: f64) -> f64 { f64::exp2(x) }
pub extern fn _log2(x: f64) -> f64 { f64::log2(x) }

fn exp2(x: Value, /*Constants{_m126_99999, _1_2, _127, exp, ..}*/_: &Constants, f: &mut FunctionBuilder<'t>) -> Value {
	/*use types::F32;
	let x = f![f fdemote(F32, x)];
	let x = f![f fmax(x, *_m126_99999)];
	let ipart = f![f fcvt_to_sint(I32, f![f fsub(x, *_1_2)])];
	let fpart = f![f fsub(x, f![f fcvt_from_sint(F32, ipart)])];
	let expipart = f![f bitcast(F32, f![f ishl_imm(f![f iadd(ipart, *_127)], 23)])];
	//let expfpart = fma![f (fma![f (exp[2], fpart, exp[1])], fpart, exp[0])];
	let expfpart = fma![f (fma![f (fma![f (fma![f (fma![f (exp[5], fpart, exp[4])], fpart, exp[3])], fpart, exp[2])], fpart, exp[1])], fpart, exp[0])];
	f![f fpromote(F64, f![f fmul(expipart, expfpart)])]*/
	let call = f![f call_indirect(
		f.import_signature(cranelift_codegen::ir::Signature{params: vec![AbiParam::new(F64)], returns: vec![AbiParam::new(F64)], /*calling_convention*/call_conv: cranelift_codegen::isa::CallConv::SystemV}),
		f![f iconst(types::I64, _exp2 as *const fn(f64)->f64 as i64)],
		&[x])];
	f.func.dfg.first_result(call)
}

fn log2(x: Value, /*Constants{_1_f32: _1, exponent, _127, mantissa, log, ..}*/_: &Constants, f: &mut FunctionBuilder<'t>) -> Value {
	/*use types::F32;
	let x = f![f fdemote(F32, x)];
	let i = f![f bitcast(I32, x)];
	let e = f![f fcvt_from_sint(F32, f![f isub(f![f ushr_imm(f![f band(i, *exponent)], 23)], *_127)])];
	let m = f![f bor(f![f bitcast(F32, f![f band(i, *mantissa)])], *_1)];
	let p = fma![f (fma![f (fma![f (fma![f (fma![f (log[5], m, log[4])], m, log[3])], m, log[2])], m, log[1])], m, log[0])];
	let p = f![f fmul(p, f![f fsub(m, *_1)])]; //?
	f![f fpromote(F64, f![f fadd(p, e)])]*/
	let call = f![f call_indirect(
		f.import_signature(cranelift_codegen::ir::Signature{params: vec![AbiParam::new(F64)], returns: vec![AbiParam::new(F64)], /*calling_convention*/call_conv: cranelift_codegen::isa::CallConv::SystemV}),
		f![f iconst(types::I64, _log2 as *const fn(f64)->f64 as i64)],
		&[x])];
	f.func.dfg.first_result(call)
}

use std::f64::consts::LN_2;

struct T { log: Value, rcp: Value, _1: Value, _2: Value, _4: Value, m: Value, mrcp: Value, rcp2: Value }

// A.T^β.exp(-Ea/kT) = exp(-(Ea/k)/T+β.logT+logA) = exp2(-(Ea/k)/T+β.log2T+log2A)
fn arrhenius(RateConstant{preexponential_factor, temperature_exponent, activation_temperature}: RateConstant, T: &T, C: &mut Constants, f: &mut FunctionBuilder) -> Value {
	if [0.,-1.,1.,2.,4.,-2.].contains(&temperature_exponent) && activation_temperature == 0. {
		let A = C.c(f, preexponential_factor);
		if temperature_exponent == 0. { A }
		else if temperature_exponent == -1. { f![f fmul(A, T.rcp)] }
		else if temperature_exponent == 1. { f![f fmul(A, T._1)] }
		else if temperature_exponent == 2. { f![f fmul(A, T._2)] }
		else if temperature_exponent == 4. { f![f fmul(A, T._4)] }
		else if temperature_exponent == -2. { f![f fmul(A, T.rcp2)] }
		else { unreachable!() }
	} else {
		let logA = C.c(f, f64::log2(preexponential_factor));
		let βlogT𐊛logA = if temperature_exponent == 0. { logA } else { fma![f (C.c(f, temperature_exponent), T.log, logA)] };
		let log_arrhenius = if activation_temperature == 0. { βlogT𐊛logA } else { fma![f (C.c(f, -activation_temperature/LN_2), T.rcp, βlogT𐊛logA)] };
		exp2(log_arrhenius, C, f)
	}
}

fn fdot<'t>(iter: impl IntoIterator<Item=(Value, impl FnOnce(&mut Constants, &mut FunctionBuilder)->Value)>, C: &mut Constants, f: &mut FunctionBuilder) -> Value {
	let mut iter = iter.into_iter();
	let mut sum = {let (a,b) = iter.next().unwrap(); f![f fmul(a, b(C, f))]};
	for (a,b) in iter { sum = fma![f (a, b(C, f), sum)]; }
	sum
}

fn dot<T>(iter: impl IntoIterator<Item=(T, Value)>, mut sum: Option<Value>, C: &mut Constants, f: &mut FunctionBuilder) -> Option<Value>
where T: num::IsZero + num::IsOne + num::IsMinusOne + Into<f64> {
	for (c,v) in iter.into_iter() {
		if c.is_zero() {}
		else if c.is_one() { sum = Some(match sum { Some(sum) => f![f fadd(sum, v)], None => v}); }
		else if c.is_minus_one() { sum = Some(match sum { Some(sum) => f![f fsub(sum, v)], None => f![f fneg(v)]}); } // fixme: reorder -a+b -> b-a to avoid neg
		else { let c = C.c(f, c.into()); sum = Some(match sum { Some(sum) => fma![f (c, v, sum)], None => f![f fmul(c, v)] }); }
	}
	sum
}

// <=> exp(dot(log(a), log(b)))
/*fn product_of_exponentiations<T>(iter: impl IntoIterator<Item=(T, Value)>, mut product: Option<Value>, C: &mut Constants, f: &mut FunctionBuilder<'t>) -> Option<Value>
where T: Into<i16> {
	for (c,v) in iter.into_iter() {
		let c = c.into();
		if c > 0 { for _ in 0..c { product = Some(match product { Some(product) => f![f fmul(product, v)], None => v}); } }
		else if c < 0 {
			let mut term = v;
			for _ in 1..-c { term = f![f fmul(term, v)]; }
			product = Some(match product { Some(product) => f![f fdiv(product, term)], None => f![f fdiv(C._1, term)]}); // fixme: reorder 1/a*b -> b/a to avoid rcp
		}
		else { assert!(c==0); }
	}
	product
}*/
fn product_of_exponentiations<T: Into<i16>>(iter: impl IntoIterator<Item=(T, Value)>, C: &mut Constants, f: &mut FunctionBuilder<'t>) -> Option<Value> {
	let (num, div) = iter.into_iter().map(|(c,v)| (c.into(), v)).filter(|&(c,_)| c!=0).partition(|&(c,_)| c>0):(Vec::<_>,Vec::<_>);
	let num = num.into_iter().fold(None, |mut a, (c,v)|{ for _ in 0..c { a = Some(match a { Some(a) => f![f fmul(a, v)], None => v }); } a });
	let div = div.into_iter().fold(None, |mut a, (c,v)|{ for _ in 0..-c { a = Some(match a { Some(a) => f![f fmul(a, v)], None => v }); } a });
	match (num, div) {
		(None, None) => None,
		(Some(num), None) => Some(num),
		(None, Some(div)) => Some(f![f fdiv(C._1, div)]), // Perf: rcp would be faster (12bit)
		(Some(num), Some(div)) => Some(f![f fdiv(num, div)]) // Perf: mul rcp would be faster (12bit)
	}
}

impl ReactionModel {
fn efficiency(&self, f: &mut FunctionBuilder<'t>, C: &mut Constants, T: &T, concentrations: &[Value], k_inf: Value) -> Value {
	use ReactionModel::*; match self {
		Elementary|Irreversible => C._1,
		ThreeBody{efficiencies} => { dot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None, C, f).unwrap() },
		PressureModification{efficiencies, k0} => {
			let Pr = f![f fmul(dot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None, C, f).unwrap(), f![f fdiv(arrhenius(*k0, T, C, f), k_inf)])];
			f![f fdiv(Pr, f![f fadd(C._1, Pr)])]
		}
		Falloff{efficiencies, k0, troe} => {
			let Pr = f![f fmul(dot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None, C, f).unwrap(), f![f fdiv(arrhenius(*k0, T, C, f), k_inf)])];
			let model::Troe{A, T3, T1, T2} = *troe;
			fn rcp(x: f64) -> f64 { 1./x }
			let Fcent = fma![f (C.c(f, 1.-A), exp2(f![f fmul(T.m, C.c(f, rcp(LN_2*T3)))], C, f),
														 fma![f (C.c(f, A), exp2(f![f fmul(T.m, C.c(f, rcp(LN_2*T1)))], C, f),
																																							exp2(f![f fmul(T.mrcp, C.c(f, T2/LN_2))], C, f) )] )];
			let logFcent = log2(Fcent, C, f);
			let c =fma![f (C.c(f, -0.67), logFcent, C.c(f, -0.4*f64::log2(10.)))];
			let N = fma![f (C.c(f, -1.27), logFcent, C.c(f, 0.75*f64::log2(10.)))];
			let logPr𐊛c = f![f fadd(log2(Pr, C, f), c)];
			let f1 = f![f fdiv(logPr𐊛c, fma![f (C.c(f, -0.14), logPr𐊛c, N)])];
			let F = exp2(f![f fdiv(logFcent, fma![f (f1, f1, C._1)])], C, f);
			f![f fmul(f![f fdiv(Pr, f![f fadd(C._1, Pr)])], F)]
		}
	}
}
}

use {binemit::CodeOffset, ir::SourceLoc};
pub struct Trap {
	pub code_offset: CodeOffset,
	pub source_location: SourceLoc,
	pub trap_code: TrapCode,
}

pub trait Rate<const CONSTANT: Property> = Fn(Constant<CONSTANT>, &StateVector<CONSTANT>, &mut Derivative<CONSTANT>, &mut [f64]);
impl Model {
pub fn rate<const CONSTANT: Property>(&self) -> (extern fn(f64, *const f64, *mut f64, *mut f64), impl Rate<CONSTANT>) {
	let mut module = cranelift_jit::JITModule::new({
		//cranelift_jit::JITBuilder::new(cranelift_module::default_libcall_names()))
		let mut flag_builder = settings::builder();
		//flag_builder.set("use_colocated_libcalls", "false").unwrap();
		//flag_builder.set("is_pic", "false").unwrap(); // FIXME set back to true once the x64 backend supports it.
		flag_builder.enable("is_pic").unwrap();
		flag_builder.set("enable_probestack", "false").unwrap();
		let isa_builder = cranelift_native::builder().unwrap();
		let isa = isa_builder.finish(settings::Flags::new(flag_builder));
		cranelift_jit::JITBuilder::with_isa(isa, cranelift_module::default_libcall_names())
	});
  let mut context = module.make_context();
	let mut function_builder_context = FunctionBuilderContext::new();
	let PTR = module.target_config().pointer_type();
	context.func.signature.params = vec![AbiParam::new(F64), AbiParam::new(PTR), AbiParam::new(PTR), AbiParam::new(PTR)];
	let mut builder = FunctionBuilder::new(&mut context.func, &mut function_builder_context);
	let entry_block = builder.create_block();
	builder.append_block_params_for_function_params(entry_block);
	builder.switch_to_block(entry_block);
	builder.seal_block(entry_block);
	let [constant, state, rate, debug]: [Value; 4] = builder.block_params(entry_block).try_into().unwrap();
	let flags = MemFlags::new();
	let ref mut f = builder;
	let ref mut C = Constants::new(f);
	let _m1 = C.c(f, -1.);
	use std::mem::size_of;
	let T = f![f load(F64, flags, state, 0*size_of::<f64>() as i32)];
	let logT = log2(T, C, f);
	let rcpT = f![f fdiv(C._1, T)];
	let T2 = f![f fmul(T, T)];
	let T3 = f![f fmul(T2, T)];
	let T4 = f![f fmul(T3, T)];
	let mT = f![f fneg(T)];
	let mrcpT = f![f fneg(rcpT)];
	let rcpT2 = f![f fmul(rcpT, rcpT)];
	let Self{species: Species{molar_mass: W, thermodynamics, heat_capacity_ratio, ..}, reactions} = self;
	let len = self.len();
	let a = thermodynamics.iter().map(|s| s.0[1]).collect(): Box<_>;
	let exp_G_RT = a[..len-1].iter().map(|a| exp2(dot(
			IntoIter::new([(a[5]/LN_2, rcpT), (-a[0], logT), (-a[1]/2./LN_2, T), ((1./3.-1./2.)*a[2]/LN_2, T2), ((1./4.-1./3.)*a[3]/LN_2, T3), ((1./5.-1./4.)*a[4]/LN_2, T4)]),
			Some(C.c(f, (a[0]-a[6])/LN_2)), C, f
		).unwrap(), C, f) ).collect():Box<_>;
	let P0_RT = f![f fmul(C.c(f, NASA7::reference_pressure), rcpT)];
	let variable = f![f load(F64, flags, state, 1*size_of::<f64>() as i32)];
	let (pressure, volume) = {use Property::*; match CONSTANT {Pressure => (constant, variable), Volume => (variable, constant)}};
	let kT = f![f fmul(C.c(f, K), T)];
	let total_concentration = f![f fdiv(pressure/*/Na*/, kT)]; // n/V = P/RT
	let amounts = (0..len-1).map(|i| f![f load(F64, flags, state, ((2+i)*size_of::<f64>()) as i32)]).collect(): Box<_>;
	let rcpV = f![f fdiv(C._1, volume)];
	let concentrations = amounts.iter().map(|&n| f![f fmul(n/*.max(0.)*/, rcpV)]).collect(): Box<_>;
	let Ca = f![f fsub(total_concentration, dot(std::iter::repeat(1.).zip(concentrations.iter().copied()), None, C, f).unwrap())];
	//if Ca < 0. { dbg!(T, C, concentrations, Ca); throw!(); }
	let ref concentrations = [&concentrations as &[_],&[Ca]].concat();
	let mut dtω = (0..len-1).map(|_| None).collect(): Box<_>;
	for (_reaction_index, Reaction{reactants, products, net, Σnet, rate_constant, model, ..}) in reactions.iter().enumerate() {
		let ref T = T{log: logT, rcp: rcpT, _1: T, _2: T2, _4: T4, m: mT, mrcp: mrcpT, rcp2: rcpT2};
		let k_inf = arrhenius(*rate_constant, T, C, f);
		let c = f![f fmul(k_inf, model.efficiency(f, C, T, concentrations, k_inf))]; // todo: CSE
		let Rf = product_of_exponentiations(reactants.iter().copied().zip(concentrations.iter().copied()), C, f).unwrap();
		let R = if let ReactionModel::Irreversible = model { Rf } else {
			let rcp_equilibrium_constant = product_of_exponentiations(net.iter().chain(&[-Σnet]).copied().zip(exp_G_RT.iter().chain(&[P0_RT]).copied()), C, f).unwrap();
			f![f store(flags, rcp_equilibrium_constant, debug, ((325+_reaction_index)*size_of::<f64>()) as i32)];
			let Rr = f![f fmul(rcp_equilibrium_constant, product_of_exponentiations(products.iter().copied().zip(concentrations.iter().copied()), C, f).unwrap())];
			f![f fsub(Rf, Rr)]
		};
		let cR = f![f fmul(c, R)];
		f![f store(flags, cR, debug, (_reaction_index*size_of::<f64>()) as i32)];
		for (index, &ν) in net.iter().enumerate() {
			let dtω = &mut dtω[index];
			match ν {
					0 => {},
					1 => match dtω { None => *dtω = Some(cR), Some(dtω) => *dtω = f![f fadd(*dtω, cR)] }
					-1 => match dtω { None => *dtω = Some(f![f fneg(cR)]), Some(dtω) => *dtω = f![f fsub(*dtω, cR)] }
					ν => match dtω { None => *dtω = Some(f![f fmul(C.c(f, ν as f64), cR)]), Some(dtω) => *dtω = fma![f (C.c(f, ν as f64), cR, *dtω)] }
			}
		}
	}

	let a = {use Property::*; match CONSTANT {Pressure => a, Volume => a.iter().zip(heat_capacity_ratio.iter()).map(|(a,γ)| a.map(|a| a / γ)).collect()}};
	struct Dot<const N: usize>(f64, [(f64, Value); N]);
	impl<'t, const N: usize> FnOnce<(&mut Constants, &mut FunctionBuilder<'_>,)> for Dot<N> {
		type Output = Value;
		extern "rust-call" fn call_once(self, (C, f,): (&mut Constants, &mut FunctionBuilder,)) -> Self::Output {
			dot(IntoIter::new(self.1), Some(C.c(f, self.0)), C, f).unwrap()
		}
	}
	let Cc/*Cp|Cv*/ = a.iter().map(|a| Dot(a[0], [(a[1], T), (a[2], T2), (a[3], T3), (a[4], T4)]));
	let m_rcp_ΣCCc = f![f fdiv(_m1, fdot(concentrations.iter().copied().zip(Cc), C, f))];
	let E_RT/*H/RT|U/RT*/ = a.iter().map(|a| Dot(a[0], [(a[5], rcpT), (a[1]/2., T), (a[2]/3., T2), (a[3]/4., T3), (a[4]/5., T4)]));
	let dtω = dtω.into_iter().map(|dtω| dtω.unwrap()).collect(): Box<_>;
	let dtT_T = f![f fmul(m_rcp_ΣCCc, fdot(dtω.iter().copied().zip(E_RT), C, f))];
	f![f store(flags, f![f fmul(dtT_T, T)], rate, 0*size_of::<f64>() as i32)];
	fn rcp(f: &mut FunctionBuilder<'t>, C: &mut Constants, x: Value) -> Value { f![f fdiv(C._1, x)] }
	let R_S_Tdtn = f![f fmul(rcp(f, C, total_concentration), dot(W[0..len-1].iter().map(|w| 1. - w/W[len-1]).zip(dtω.iter().copied()), None, C, f).unwrap())]; // R/A Tdtn (constant pressure: A=V, constant volume: A=P)
	let dtS_S = f![f fadd(R_S_Tdtn, dtT_T)];
	f![f store(flags, f![f fmul(dtS_S, variable)], rate, 1*size_of::<f64>() as i32)];
	//let dtn = dtω.into_iter().map(|&dtω| f![f fmul(volume, dtω)]).collect(): Box<_>;
	//for (i, &dtn) in dtn.into_iter().enumerate() { f![f store(flags, dtn, rate, ((2+i)*size_of::<f64>()) as i32)]; }
	for (i, &dtω) in dtω.into_iter().enumerate() { f![f store(flags, f![f fmul(volume, dtω)], rate, ((2+i)*size_of::<f64>()) as i32)]; }
	builder.ins().return_(&[]);
	builder.finalize();
	let clif = builder.display(None);
	//eprintln!("{}", clif);
	//std::fs::write("/tmp/CH4+O2.clif", clif.to_string()).unwrap();
	/*if false {
		let function = cranelift_reader::parse_functions(&clif.to_string()).unwrap().remove(0);
		let mut context = codegen::Context::new();
		context.func = function;
		let mut mem = vec![];
		let isa_builder = isa::lookup(target_lexicon::Triple::host()).unwrap();
		let mut flag_builder = settings::builder();
		flag_builder.enable("is_pic").unwrap();
		flag_builder.set("enable_probestack", "false").unwrap();
		let isa = isa_builder.finish(settings::Flags::new(flag_builder));
		let code_info = context.compile_and_emit(&*isa, &mut mem, &mut binemit::NullRelocSink{}, &mut binemit::NullTrapSink{}, &mut binemit::NullStackMapSink{});
		use capstone::arch::BuildsCapstone;
		let capstone = capstone::Capstone::new().x86().mode(capstone::arch::x86::ArchMode::Mode64).build().unwrap();
		let instructions = capstone.disasm_all(&mem, 0).unwrap();
		let mut file = std::fs::File::create("/tmp/CH4+O2.asm").unwrap();
		for i in instructions.iter() { use std::io::Write; writeln!(file, "{}\t{}", i.mnemonic().unwrap(), i.op_str().unwrap()).unwrap(); }
	}*/
  let id = module.declare_function(&"", Linkage::Export, &context.func.signature).unwrap();
  module.define_function(id, &mut context, &mut binemit::NullTrapSink{}).unwrap();
	module.finalize_definitions();
	let function = module.get_finalized_function(id);
	let function = unsafe{std::mem::transmute::<_,extern fn(f64, *const f64, *mut f64, *mut f64)>(function)};
	(function, move |constant:Constant<CONSTANT>, state:&StateVector<CONSTANT>, derivative:&mut Derivative<CONSTANT>, debug: &mut [f64]| {
		let constant = constant.0;
		function(constant, state.0.as_ptr(), derivative.0.as_mut_ptr(), debug.as_mut_ptr());
	})
}
}