use std::default::default;
use super::*;
pub use cranelift::codegen::ir::{function::Function, types::{I32, I64, F32, F64}};
use cranelift::{
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

//#[derive(Debug, num_enum::TryFromPrimitive)] #[allow(non_camel_case_types)] #[repr(u32)] pub enum Intrinsic { exp2, log2 }

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
  fn load(&mut self, base: Value, index: usize) -> Value { self.ins().load(F64, MemFlags::trusted(), base, (index*std::mem::size_of::<f64>()) as i32) }
	fn store(&mut self, value: Value, base: Value, index: usize) { self.ins().store(MemFlags::trusted(), value, base, (index*std::mem::size_of::<f64>()) as i32); }
	fn neg(&mut self, x: Value) -> Value { self.ins().fneg(x) }
	fn min(&mut self, x: Value, y: Value) -> Value { self.ins().fmin(x, y) }
	fn max(&mut self, x: Value, y: Value) -> Value { self.ins().fmax(x, y) }
	fn add(&mut self, x: Value, y: Value) -> Value { self.ins().fadd(x, y) }
	fn sub(&mut self, x: Value, y: Value) -> Value { self.ins().fsub(x, y) }
	fn mul(&mut self, x: Value, y: Value) -> Value { self.ins().fmul(x, y) }
	fn div(&mut self, x: Value, y: Value) -> Value { self.ins().fdiv(x, y) }
	fn fma(&mut self, x: Value, y: Value, z: Value) -> Value { let mul = self.mul(x, y); self.add(mul, z) }
}

// Evaluate arguments before borrowing builder for main instruction
fn store(value: Value, base: Value, index: usize, f: &mut FunctionBuilder<'_>) { f.store(value, base, index); }
fn max(x: Value, y: Value, f: &mut FunctionBuilder<'_>) -> Value { f.max(x, y) }
fn add(x: Value, y: Value, f: &mut FunctionBuilder<'_>) -> Value { f.add(x, y) }
fn sub(x: Value, y: Value, f: &mut FunctionBuilder<'_>) -> Value { f.sub(x, y) }
fn mul(x: Value, y: Value, f: &mut FunctionBuilder<'_>) -> Value { f.mul(x, y) }
fn div(x: Value, y: Value, f: &mut FunctionBuilder<'_>) -> Value { f.div(x, y) }
fn fma(x: Value, y: Value, z: Value, f: &mut FunctionBuilder<'_>) -> Value { f.fma(x, y, z) }

impl FunctionBuilder<'_> {
	fn rcp(&mut self, x: Value) -> Value { div(self.f64(1.), x, self) }
	fn c(&mut self, value: f64) -> Value { self.f64(value) }
	fn dot(&mut self, iter: impl IntoIterator<Item=(Value, Value)>) -> Value {
		let mut iter = iter.into_iter();
		let mut sum = { let (a,b) = iter.next().unwrap(); self.mul(a, b) };
		for (a,b) in iter { sum = self.fma(a, b, sum); }
		sum
	}
	fn fdot(&mut self, iter: impl IntoIterator<Item=(Value, impl FnOnce(&mut FunctionBuilder)->Value)>) -> Value {
		let mut iter = iter.into_iter();
		let mut sum = { let (a,b) = iter.next().unwrap(); let b = b(self); self.mul(a, b) };
		for (a,b) in iter { let b = b(self); sum = self.fma(a, b, sum); }
		sum
	}
	fn cdot<T>(&mut self, iter: impl IntoIterator<Item=(T, Value)>, mut sum: Option<Value>) -> Option<Value>
	where T: num::IsZero + num::IsOne + num::IsMinusOne + Into<f64> {
		let f = self;
		for (c,v) in iter.into_iter() {
			if c.is_zero() {}
			else if c.is_one() { sum = Some(match sum { Some(sum) => f.add(sum, v), None => v}); }
			else if c.is_minus_one() { sum = Some(match sum { Some(sum) => f.sub(sum, v), None => f.neg(v)}); } // fixme: reorder -a+b -> b-a to avoid neg
			else { let c = f.c(c.into()); sum = Some(match sum { Some(sum) => f.fma(c, v, sum), None => f.mul(c, v) }); }
		}
		sum
	}
	fn product_of_exponentiations<T: Into<i16>>(&mut self, iter: impl IntoIterator<Item=(T, Value)>) -> Option<Value> {
		let f = self;
		let (num, div) : (Vec::<_>,Vec::<_>) = iter.into_iter().map(|(c,v)| (c.into(), v)).filter(|&(c,_)| c!=0).partition(|&(c,_)| c>0);
		let num = num.into_iter().fold(None, |mut a, (c,v)|{ for _ in 0..c { a = Some(match a { Some(a) => f.mul(a, v), None => v }); } a });
		let div = div.into_iter().fold(None, |mut a, (c,v)|{ for _ in 0..-c { a = Some(match a { Some(a) => f.mul(a, v), None => v }); } a });
		match (num, div) {
			(None, None) => None,
			(Some(num), None) => Some(num),
			(None, Some(div)) => Some(f.rcp(div)),
			(Some(num), Some(div)) => Some(f.div(num, div)) // Perf: mul rcp would be faster (12bit)
		}
	}
}

macro_rules! f { // Evaluates args before borrowing self
	[$f:ident $function:ident($arg0:expr)] => {{
		let arg0 = $arg0;
		$f.ins().$function(arg0)
	}};
	[$f:ident $function:ident($arg0:expr, $arg1:expr)] => {{
		let arg0 = $arg0;
		let arg1 = $arg1;
		$f.ins().$function(arg0, arg1)
	}};
}

impl FunctionBuilder<'_> {
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
}

fn exp2(x: Value, f: &mut FunctionBuilder<'_>) -> Value { f.exp2(x) }

use std::f64::consts::LN_2;

struct T { log: Value, rcp: Value, _1: Value, _2: Value, _4: Value, m: Value, mrcp: Value, rcp2: Value }

// A.T^β.exp(-Ea/kT) = exp(-(Ea/k)/T+β.logT+logA) = exp2(-(Ea/k)/T+β.log2T+log2A)
fn arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant, T: &T, f: &mut FunctionBuilder) -> Value {
	if [0.,-1.,1.,2.,4.,-2.].contains(&temperature_exponent) && activation_temperature == 0. {
		let A = f.c(A);
		if temperature_exponent == 0. { A }
		else if temperature_exponent == -1. { f.mul(A, T.rcp) }
		else if temperature_exponent == 1. { f.mul(A, T._1) }
		else if temperature_exponent == 2. { f.mul(A, T._2) }
		else if temperature_exponent == 4. { f.mul(A, T._4) }
		else if temperature_exponent == -2. { f.mul(A, T.rcp2) }
		else { unreachable!() }
	} else {
		let βlogT𐊛logA = if temperature_exponent == 0. { f.c(f64::log2(A)) } else { fma(f.c(temperature_exponent), T.log, f.c(f64::log2(A)), f) };
		let log_arrhenius = if activation_temperature == 0. { βlogT𐊛logA } else { fma(f.c(-activation_temperature/LN_2), T.rcp, βlogT𐊛logA, f) };
		f.exp2(log_arrhenius)
	}
}

impl ReactionModel {
fn efficiency(&self, T: &T, concentrations: &[Value], k_inf: Value, f: &mut FunctionBuilder) -> Value {
	use ReactionModel::*; match self {
		Elementary|Irreversible => f.c(1.),
		ThreeBody{efficiencies} => f.cdot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None).unwrap(),
		PressureModification{efficiencies, k0} => {
			let Pr = mul(f.cdot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None).unwrap(), div(arrhenius(k0, T, f), k_inf, f), f);
			div(Pr, add(f.c(1.), Pr, f), f)
		}
		Falloff{efficiencies, k0, troe} => {
			let Pr = mul(f.cdot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None).unwrap(), div(arrhenius(k0, T, f), k_inf, f), f);
			let model::Troe{A, T3, T1, T2} = *troe;
			let Fcent = fma(f.c(1.-A), exp2(mul(T._1, f.c(-1./(LN_2*T3)), f), f),
												 fma(f.c(A), exp2(mul(T._1, f.c(-1./(LN_2*T1)), f), f),
												                        exp2(mul(T.rcp, f.c(-T2/LN_2), f), f), f), f);
			let logFcent = f.log2(Fcent);
			let c =fma(f.c(-0.67), logFcent, f.c(-0.4*f64::log2(10.)), f);
			let N = fma(f.c(-1.27), logFcent, f.c(0.75*f64::log2(10.)), f);
			let logPr𐊛c = add(f.log2(Pr), c, f);
			let f1 = div(logPr𐊛c, fma(f.c(-0.14), logPr𐊛c, N, f), f);
			let F = exp2(div(logFcent, fma(f1, f1, f.c(1.), f), f), f);
			mul(div(Pr, add(f.c(1.), Pr, f), f), F, f)
		}
	}
}
}

fn dot<const N: usize>(constant: f64, iter: [(f64, Value); N]) -> impl /*FnOnce<(&'t mut FunctionBuilder<'t>,)>*/FnOnce(&mut FunctionBuilder<'_>)->Value {
	move |f: &mut FunctionBuilder<'_>| { let c = f.c(constant); f.cdot(IntoIter::new(iter), Some(c)).unwrap() }
}

fn reaction<'t, Reactions: IntoIterator<Item=&'t Reaction>>(Species{thermodynamics, ..}: &Species, reactions: Reactions, concentrations: &[Value], f: &mut FunctionBuilder, T: Value, rcp_T: Value, T2: Value, T3: Value, T4: Value) -> Box<[Value]> {
	let logT = f.log2(T);
	let mT = f.neg(T);
	let mrcp_T = f.neg(rcp_T);
	let rcp_T2 = f.mul(rcp_T, rcp_T);
	let a = map(&**thermodynamics, |s| s.pieces[1]);
	let exp_Gibbs_RT = map(&a[..a.len()-1], |a|
		exp2(dot((a[0]-a[6])/LN_2, [(a[5]/LN_2, rcp_T), (-a[0], logT), (-a[1]/2./LN_2, T), ((1./3.-1./2.)*a[2]/LN_2, T2), ((1./4.-1./3.)*a[3]/LN_2, T3), ((1./5.-1./4.)*a[4]/LN_2, T4)])(f), f)
	);
	let P0_RT = mul(f.c(NASA7::reference_pressure), rcp_T, f); // TODO: CSE P0_RT^-Σnet
	let mut dtω = vec![None; a.len()-1].into_boxed_slice();
	for (_reaction_index, reaction) in reactions.into_iter().enumerate() {
		let Reaction{reactants, products, net, Σnet, rate_constant, model, ..} = reaction;
		let ref T = T{log: logT, rcp: rcp_T, _1: T, _2: T2, _4: T4, m: mT, mrcp: mrcp_T, rcp2: rcp_T2};
		let k_inf = arrhenius(rate_constant, T, f);
		let c = mul(k_inf, model.efficiency(T, &concentrations, k_inf, f), f); // todo: CSE
		let Rf = f.product_of_exponentiations(reactants.iter().copied().zip(concentrations.iter().copied())).unwrap();
		let R = if let ReactionModel::Irreversible = model { Rf } else {
			let rcp_equilibrium_constant = f.product_of_exponentiations(net.iter().chain(&[-Σnet]).copied().zip(exp_Gibbs_RT.iter().chain(&[P0_RT]).copied())).unwrap();
			let Rr = mul(rcp_equilibrium_constant, f.product_of_exponentiations(products.iter().copied().zip(concentrations.iter().copied())).unwrap(), f);
			f.sub(Rf, Rr)
		};
		let cR = f.mul(c, R);
		for (index, &ν) in net.iter().enumerate() {
			let dtω = &mut dtω[index];
			match ν {
					0 => {},
					1 => match dtω { None => *dtω = Some(cR), Some(dtω) => *dtω = f.add(*dtω, cR) }
					-1 => match dtω { None => *dtω = Some(f.neg(cR)), Some(dtω) => *dtω = f.sub(*dtω, cR) }
					ν => match dtω { None => *dtω = Some(mul(f.c(ν as f64), cR, f)), Some(dtω) => *dtω = fma(f.c(ν as f64), cR, *dtω, f) }
			}
		}
	}
	map(&*dtω, |dtω| dtω.unwrap())
}

pub fn rate<'t, Reactions: IntoIterator<Item=&'t Reaction>, const CONSTANT: Property>(species@Species{molar_mass, thermodynamics, heat_capacity_ratio, ..}: &Species, reactions: Reactions, stride: usize) -> Function {
	let mut function = Function::new();
	let mut function_builder_context = FunctionBuilderContext::new();
	let ref mut f = FunctionBuilder::new(&mut function, &mut function_builder_context);
	let entry_block = f.create_block();
	const parameters: [Type; 4] = [I64, I64, I64, F64];
	f.func.signature.params = map(&parameters, |t| AbiParam::new(*t)).to_vec();
	f.append_block_params_for_function_params(entry_block);
	f.switch_to_block(entry_block);
	f.seal_block(entry_block);
	let [index, states, rates, constant]: [Value; 4] = f.block_params(entry_block).try_into().unwrap();
	let offset = f.ins().ishl_imm(index, 3);
	let state = f.ins().iadd(states, offset);
	let rates = f.ins().iadd(rates, offset);
	let T = f.load(states, 0*stride);
	let rcp_T = f.rcp(T);
	let variable = f.load(state, 1*stride);
	let (pressure_R, volume) = {use Property::*; match CONSTANT {Pressure => (constant, variable), Volume => (variable, constant)}};
	let total_concentration = f.mul(pressure_R, rcp_T);
	let active_amounts = iter::eval(species.len()-1, |i| f.load(state, (2+i)*stride));
	let rcpV = f.rcp(volume);
	let active_amounts = map(&*active_amounts, |&n| max(f.c(0.), n, f));
	let active_concentrations = map(&*active_amounts, |&n| f.mul(n, rcpV));
	let inert_concentration = sub(total_concentration, f.cdot(std::iter::repeat(1.).zip(active_concentrations.iter().copied()), None).unwrap(), f);
	let ref concentrations = [&active_concentrations as &[_],&[inert_concentration]].concat();
	let T2 = f.mul(T, T);
	let T3 = f.mul(T2, T);
	let T4 = f.mul(T3, T);
	let dtω = reaction(species, reactions, &*concentrations, f, T, rcp_T, T2, T3, T4);
	let a = map(&**thermodynamics, |s| *s.piece(/*T:*/1000.));
	let a = {use Property::*; match CONSTANT {Pressure => a, Volume => a.iter().zip(heat_capacity_ratio.iter()).map(|(a,γ)| a.map(|a| a / γ)).collect()}};
	let E_RT/*H/RT|U/RT*/ = a.iter().map(|a| dot(a[0], [(a[5], rcp_T), (a[1]/2., T), (a[2]/3., T2), (a[3]/4., T3), (a[4]/5., T4)]));
	let E_RT = dtω.into_iter().enumerate().zip(E_RT).zip(&**molar_mass).map(|(((i, &dtω), E_RT), molar_mass)| move |f: &mut FunctionBuilder<'_>| {
		store(dtω, rates, (2+i)*stride, f);
		E_RT(f)
	});
	let Cc/*Cp|Cv*/ = a.iter().map(|a| dot(a[0], [(a[1], T), (a[2], T2), (a[3], T3), (a[4], T4)]));
	let m_rcp_ΣCCc = div(f.c(-1.), f.fdot(concentrations.iter().copied().zip(Cc)), f);
	macro_rules! dot { ($a:ident, $b:ident, $f:ident) => ($f.fdot($a.iter().copied().zip($b))) }
	let dtT_T = mul(m_rcp_ΣCCc, dot!(dtω, E_RT, f), f);
	store(mul(dtT_T, T, f), rates, 0*stride, f);
	let len = species.len();
	let R_S_Tdtn = mul(f.rcp(total_concentration), f.cdot(molar_mass[0..len-1].iter().map(|w| 1. - w/molar_mass[len-1]).zip(dtω.iter().copied()), None).unwrap(), f);
	let dtS_S = f.add(R_S_Tdtn, dtT_T);
	store(f.mul(dtS_S, variable), rates, 1*stride, f);
	f.ins().return_(&[]);
	function
}

pub trait Rate<const CONSTANT: Property> = Fn(Constant<CONSTANT>, &StateVector<CONSTANT>, &mut Derivative<CONSTANT>);
#[cfg(feature="jit")] pub fn JIT<'t, const CONSTANT: Property>(species: &Species, reactions: impl IntoIterator<Item=&'t Reaction>) -> impl Rate<CONSTANT> {
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
	context.func = rate::<_, CONSTANT>(species, reactions, 1);
  let id = module.declare_function(&"", cranelift_module::Linkage::Export, &context.func.signature).unwrap();
  module.define_function(id, &mut context, &mut cranelift_codegen::binemit::NullTrapSink{}, &mut cranelift_codegen::binemit::NullStackMapSink{}).unwrap();
	module.finalize_definitions();
	let function = module.get_finalized_function(id);
	let function = unsafe{std::mem::transmute::<_,extern fn(u64, f64, *const f64, *mut f64)>(function)};
	move |constant:Constant<CONSTANT>, state:&StateVector<CONSTANT>, derivative:&mut Derivative<CONSTANT>| {
		let constant = constant.0;
		function(0, constant, state.0.as_ptr(), derivative.0.as_mut_ptr());
	}
}
