use std::{default::default, convert::TryFrom};
use super::*;

#[derive(Clone, Copy)] pub struct RateConstant {
	pub preexponential_factor: f64,
	pub temperature_exponent: f64,
	pub activation_temperature: f64
}

impl From<&model::RateConstant> for RateConstant {
	fn from(model::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: &model::RateConstant) -> Self {
		const J_per_cal: f64 = 4.184;
		Self{preexponential_factor: *preexponential_factor, temperature_exponent: *temperature_exponent, activation_temperature: activation_energy*J_per_cal/(K*NA)}
	}
}

use model::Troe;

pub enum ReactionModel {
	Elementary,
	Irreversible,
	ThreeBody { efficiencies: Box<[f64]> },
	PressureModification { efficiencies: Box<[f64]>, k0: RateConstant },
	Falloff { efficiencies: Box<[f64]>, k0: RateConstant, troe: Troe },
}

pub struct Reaction {
	pub reactants: Box<[u8]>,
	pub products: Box<[u8]>,
	pub net: Box<[i8/*; S-1*/]>,
	pub Σreactants: u8,
	pub Σproducts: u8,
	pub Σnet: i8,
	pub rate_constant: RateConstant,
	pub model: ReactionModel,
}

use iter::{map, eval};

impl Reaction {
	pub fn new(species_names: &[&str], model::Reaction{ref equation, rate_constant, model}: &model::Reaction) -> Self {
		for side in equation { for (specie, _) in side { assert!(species_names.contains(&specie), "{}", specie) } }
		let [reactants, products] = iter::vec::eval(equation, |e| species_names.iter().map(|&s| *e.get(s).unwrap_or(&0)).collect::<Box<_>>());
		let net = products.into_iter().zip(reactants.into_iter()).take(species_names.len()-1).map(|(&a, &b)| a as i8 - b as i8).collect();
		/*{let mut net_composition = Map::new();
			for (s, &ν) in net.into_iter().enumerate() {
				for (element, &count) in species_composition[s] {
					if !net_composition.contains_key(&element) { net_composition.insert(element, 0); }
					*net_composition.get_mut(&element).unwrap() += ν as i8 * count as i8;
				}
			}
			for (_, &ν) in &net_composition { assert!(ν == 0, "{:?} {:?}", net_composition, equation); }
		}*/
		let [Σreactants, Σproducts] = [reactants.iter().sum(), products.iter().sum()];
		let Σnet = Σproducts as i8 - Σreactants as i8;
		let from = |efficiencies:&Map<_,_>| map(species_names, |&specie| *efficiencies.get(specie).unwrap_or(&1.));
		Reaction{
			reactants, products, net, Σreactants, Σproducts, Σnet,
			rate_constant: rate_constant.into(),
			model: {use model::ReactionModel::*; match model {
				Elementary => ReactionModel::Elementary,
				Irreversible => ReactionModel::Irreversible,
				ThreeBody{efficiencies} => ReactionModel::ThreeBody{efficiencies: from(efficiencies)},
				PressureModification{efficiencies, k0} => ReactionModel::PressureModification{efficiencies: from(efficiencies), k0: k0.into()},
				Falloff{efficiencies, k0, troe} => ReactionModel::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe: *troe},
			}}
		}
	}
}

#[derive(PartialEq, Eq)] pub enum Property { Pressure, Volume }
#[derive(Clone, Copy)] pub struct Constant<const CONSTANT: Property>(pub f64);
#[derive(derive_more::Deref, Default)] pub struct StateVector<const CONSTANT: Property>(pub Box<[f64/*; T,P|V,[S-1]*/]>);
pub type Derivative<const CONSTANT: Property> = StateVector<CONSTANT>;

impl State {
	//pub fn constant<const CONSTANT: Property>(&Self{pressure, volume, ..}: &Self) -> f64 { // arbitrary_self_types
	pub fn constant<const CONSTANT: Property>(&self) -> Constant<CONSTANT> { let Self{pressure_R, volume, ..} = self;
		Constant(*{use Property::*; match CONSTANT {Pressure => pressure_R, Volume => volume}})
	}
}

use std::array::IntoIter;

impl<const CONSTANT: Property> From<&State> for StateVector<CONSTANT> {
	fn from(State{temperature, pressure_R, volume, amounts}: &State) -> Self {
		Self([*temperature, *{use Property::*; match CONSTANT {Pressure => volume, Volume => pressure_R}}].iter().chain(amounts[..amounts.len()-1].iter()).copied().collect())
	}
}

impl State {
	#[track_caller] pub fn new<const CONSTANT: Property>(total_amount: f64, Constant(thermodynamic_state_constant): Constant<CONSTANT>, u: &StateVector<CONSTANT>) -> Self {
		let u = &u.0;
		//assert!(!u.iter().any(|&n| n<0.));
		let amounts = &u[2..];
		let (pressure_R, volume) = {use Property::*; match CONSTANT {
			Pressure => (thermodynamic_state_constant, u[1]),
			Volume => (u[1], thermodynamic_state_constant)
		}};
		State{
			temperature: u[0],
			pressure_R,
			volume,
			amounts: amounts.iter().chain(&[total_amount - iter::into::Sum::<f64>::sum(amounts)]).copied().collect()
		}
	}
}

impl std::fmt::Display for State {
	fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
		let Self{temperature, pressure_R, volume, amounts} = self;
		write!(fmt, "T: {}, P: {}, V: {}, n: {:?}", temperature/*/K*/, pressure_R*(K*NA), volume, amounts)
	}
}

use cranelift::{
	frontend::{self as frontend, /*FunctionBuilder,*/ FunctionBuilderContext},
	codegen::{
		ir::{function::Function, Signature, InstBuilder, types::{I64, F64}, MemFlags, entities::{Value, SigRef, FuncRef}, AbiParam, ExternalName, ExtFuncData},
		isa::CallConv,
		//write::write_function
	},
};

fn import(function: &mut Function, index: u32, signature: SigRef) {
	assert!(function.import_function(ExtFuncData{name: ExternalName::User{namespace: 0, index}, signature, colocated: default()}).as_u32() == index)
}

#[derive(derive_more::Deref,derive_more::DerefMut)] struct FunctionBuilder<'t> {
	#[deref]#[deref_mut] builder: frontend::FunctionBuilder<'t>,
	constants: std::collections::HashMap<u64, Value>,
}

#[derive(Debug, num_enum::TryFromPrimitive)] #[allow(non_camel_case_types)] #[repr(u32)] enum Intrinsic { exp2, log2 }

impl FunctionBuilder<'t> {
	fn new(function: &'t mut Function, function_builder_context: &'t mut FunctionBuilderContext) -> Self { Self{
			builder: frontend::FunctionBuilder::new(function, function_builder_context),
			constants: default(),
	} }
	fn f64(&mut self, value: f64) -> Value {
		match self.constants.entry(value.to_bits()) {
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
	fn c(&mut self, value: f64) -> Value { self.f64(value) }
	fn dot(&mut self, iter: impl IntoIterator<Item=(Value, impl FnOnce(&mut FunctionBuilder)->Value)>) -> Value {
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
}

struct Constants {
	//_0: Value,
	_1: Value,
}
impl Constants {
	fn new(f: &mut FunctionBuilder<'_>) -> Self { Self {
		//_0: f.c(0.),
		_1: f.c(1.),
	} }
}

#[derive(derive_more::Deref,derive_more::DerefMut)] struct Builder<'t> { // with intrinsics and static constants
	#[deref]#[deref_mut] builder: FunctionBuilder<'t>,
	c: Constants,
}

impl Builder<'t> {
	fn new(mut builder: FunctionBuilder<'t>) -> Self {
		let signature = builder.func.import_signature(Signature{params: vec![AbiParam::new(F64)], returns: vec![AbiParam::new(F64)], call_conv: CallConv::Fast});
		import(builder.func, Intrinsic::exp2 as u32, signature);
		import(builder.func, Intrinsic::log2 as u32, signature);
		let constants = Constants::new(&mut builder);
		Self{builder, c: constants}
	}
	fn rcp(&mut self, x: Value) -> Value { self.builder.div(self.c._1, x) }
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
	fn exp2(&mut self, x: Value) -> Value {
		let call = self.builder.ins().call(FuncRef::from_u32(Intrinsic::exp2 as u32), &[x]);
		self.func.dfg.first_result(call)
	}
	fn log2(&mut self, x: Value) -> Value {
		let call = self.builder.ins().call(FuncRef::from_u32(Intrinsic::log2 as u32), &[x]);
		self.func.dfg.first_result(call)
	}
}

fn exp2(x: Value, f: &mut Builder<'_>) -> Value { f.exp2(x) }

use std::f64::consts::LN_2;

struct T { log: Value, rcp: Value, _1: Value, _2: Value, _4: Value, m: Value, mrcp: Value, rcp2: Value }

// A.T^β.exp(-Ea/kT) = exp(-(Ea/k)/T+β.logT+logA) = exp2(-(Ea/k)/T+β.log2T+log2A)
fn arrhenius(&RateConstant{preexponential_factor, temperature_exponent, activation_temperature}: &RateConstant, T: &T, f: &mut Builder) -> Value {
	if [0.,-1.,1.,2.,4.,-2.].contains(&temperature_exponent) && activation_temperature == 0. {
		let A = f.c(preexponential_factor);
		if temperature_exponent == 0. { A }
		else if temperature_exponent == -1. { f.mul(A, T.rcp) }
		else if temperature_exponent == 1. { f.mul(A, T._1) }
		else if temperature_exponent == 2. { f.mul(A, T._2) }
		else if temperature_exponent == 4. { f.mul(A, T._4) }
		else if temperature_exponent == -2. { f.mul(A, T.rcp2) }
		else { unreachable!() }
	} else {
		let logA = f.c(f64::log2(preexponential_factor));
		let βlogT𐊛logA = if temperature_exponent == 0. { logA } else { fma(f.c(temperature_exponent), T.log, logA, f) };
		let log_arrhenius = if activation_temperature == 0. { βlogT𐊛logA } else { fma(f.c(-activation_temperature/LN_2), T.rcp, βlogT𐊛logA, f) };
		f.exp2(log_arrhenius)
	}
}

impl ReactionModel {
fn efficiency(&self, T: &T, concentrations: &[Value], k_inf: Value, f: &mut Builder) -> Value {
	use ReactionModel::*; match self {
		Elementary|Irreversible => f.c(1.),
		ThreeBody{efficiencies} => { f.cdot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None).unwrap() },
		PressureModification{efficiencies, k0} => {
			let Pr = mul(f.cdot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None).unwrap(), div(arrhenius(k0, T, f), k_inf, f), f);
			div(Pr, add(f.c._1, Pr, f), f)
		}
		Falloff{efficiencies, k0, troe} => {
			let Pr = mul(f.cdot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None).unwrap(), div(arrhenius(k0, T, f), k_inf, f), f);
			let model::Troe{A, T3, T1, T2} = *troe;
			fn rcp(x: f64) -> f64 { 1./x }
			let Fcent = fma(f.c(1.-A), exp2(mul(T.m, f.c(rcp(LN_2*T3)), f), f),
												 fma(f.c(A), exp2(mul(T.m, f.c(rcp(LN_2*T1)), f), f),
												                        exp2(mul(T.mrcp, f.c(T2/LN_2), f), f), f), f);
			let logFcent = f.log2(Fcent);
			let c =fma(f.c(-0.67), logFcent, f.c(-0.4*f64::log2(10.)), f);
			let N = fma(f.c(-1.27), logFcent, f.c(0.75*f64::log2(10.)), f);
			let logPr𐊛c = add(f.log2(Pr), c, f);
			let f1 = div(logPr𐊛c, fma(f.c(-0.14), logPr𐊛c, N, f), f);
			let F = exp2(div(logFcent, fma(f1, f1, f.c._1, f), f), f);
			mul(div(Pr, add(f.c._1, Pr, f), f), F, f)
		}
	}
}
}

use std::fmt::Write;
//pub trait Rate<const CONSTANT: Property> = Fn(Constant<CONSTANT>, &StateVector<CONSTANT>, &mut Derivative<CONSTANT>, &mut [f64]);
#[fehler::throws(std::fmt::Error)] pub fn rate<'t, Reactions: IntoIterator<Item=&'t Reaction>, const CONSTANT: Property>(species@Species{molar_mass, thermodynamics, heat_capacity_ratio, ..}: &Species, reactions: Reactions, stride: usize) -> String {
	let mut function = Function::new();
	function.signature.params = vec![/*AbiParam::new(F64),*/ AbiParam::new(I64), AbiParam::new(I64)];
	let mut function_builder_context = FunctionBuilderContext::new();
	let mut f = FunctionBuilder::new(&mut function, &mut function_builder_context);
	let entry_block = f.create_block();
	f.append_block_params_for_function_params(entry_block);
	f.switch_to_block(entry_block);
	f.seal_block(entry_block);
	let [/*constant,*/ state, rates]: [Value; 2] = f.block_params(entry_block).try_into().unwrap();
	let ref mut f = Builder::new(f);
	let T = f.load(state, 0*stride);
	let logT = f.log2(T);
	let rcpT = f.rcp(T);
	let T2 = f.mul(T, T);
	let T3 = f.mul(T2, T);
	let T4 = f.mul(T3, T);
	let mT = f.neg(T);
	let mrcpT = f.neg(rcpT);
	let rcpT2 = f.mul(rcpT, rcpT);
	let len = species.len();
	let a = map(&**thermodynamics, |s| s.0[1]);

	fn dot<const N: usize>(constant: f64, iter: [(f64, Value); N]) -> impl /*FnOnce<(&'t mut FunctionBuilder<'t>,)>*/FnOnce(&mut FunctionBuilder<'_>)->Value {
		move |f: &mut FunctionBuilder<'_>| { let c = f.c(constant); f.cdot(IntoIter::new(iter), Some(c)).unwrap() }
	}
	let exp_G_RT = map(&a[..len-1], |a|
		exp2(dot((a[0]-a[6])/LN_2, [(a[5]/LN_2, rcpT), (-a[0], logT), (-a[1]/2./LN_2, T), ((1./3.-1./2.)*a[2]/LN_2, T2), ((1./4.-1./3.)*a[3]/LN_2, T3), ((1./5.-1./4.)*a[4]/LN_2, T4)])(f), f)
	);
	let P0_RT = mul(f.c(NASA7::reference_pressure), rcpT, f);
	/*let variable = f.load(state, 1*stride);
	let (pressure_R, volume) = {use Property::*; match CONSTANT {Pressure => (constant, variable), Volume => (variable, constant)}};
	let total_concentration = f.div(pressure_R, T); // n/V = P/RT
	let amounts = eval(len-1, |i| f.load(state, (2+i)*stride));
	let rcpV = f.rcp(volume);
	let amounts = map(&*amounts, |&n| max(f.c._0, n, f));
	let concentrations = map(&*amounts, |&n| f.mul(n, rcpV));
	let Ca = sub(total_concentration, f.cdot(std::iter::repeat(1.).zip(concentrations.iter().copied()), None).unwrap(), f);
	let ref concentrations = [&concentrations as &[_],&[Ca]].concat();*/
	let concentrations = eval(len-1, |i| f.load(state, (1+i)*stride));
	let mut dtω = vec![None; len-1].into_boxed_slice();
	for (_reaction_index, reaction) in reactions.into_iter().enumerate() {
		let Reaction{reactants, products, net, Σnet, rate_constant, model, ..} = reaction;
		let ref T = T{log: logT, rcp: rcpT, _1: T, _2: T2, _4: T4, m: mT, mrcp: mrcpT, rcp2: rcpT2};
		let k_inf = arrhenius(rate_constant, T, f);
		let c = mul(k_inf, model.efficiency(T, &concentrations, k_inf, f), f); // todo: CSE
		let Rf = f.product_of_exponentiations(reactants.iter().copied().zip(concentrations.iter().copied())).unwrap();
		let R = if let ReactionModel::Irreversible = model { Rf } else {
			let rcp_equilibrium_constant = f.product_of_exponentiations(net.iter().chain(&[-Σnet]).copied().zip(exp_G_RT.iter().chain(&[P0_RT]).copied())).unwrap();
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

	let a = {use Property::*; match CONSTANT {Pressure => a, Volume => a.iter().zip(heat_capacity_ratio.iter()).map(|(a,γ)| a.map(|a| a / γ)).collect()}};
	let E_RT/*H/RT|U/RT*/ = a.iter().map(|a| dot(a[0], [(a[5], rcpT), (a[1]/2., T), (a[2]/3., T2), (a[3]/4., T3), (a[4]/5., T4)]));
	let dtω = map(&*dtω, |dtω| dtω.unwrap());
	//for (i, &dtω) in dtω.into_iter().enumerate() { store(f.mul(volume, dtω), rates, (/*2*/+i)*stride, f); }
	let E_RT = dtω.into_iter().enumerate().zip(E_RT).map(|((i, &dtω), E_RT)| move |f: &mut FunctionBuilder<'_>| {
		store(dtω, rates, (/*2*/1+i)*stride, f); // Store dtω when loading dtω for dtω.E_RT
		E_RT(f)
	});
	//fn dot(array: &[Value], iter: impl Iterator<Item=FnOnce(&mut FunctionBuilder)->Value>, f: &mut FunctionBuilder) -> Value { f.dot(a.iter().copied().zip(b)) }
	macro_rules! dot { ($a:ident, $b:ident, $f:ident) => ($f.dot($a.iter().copied().zip($b))) }
	let heat_release_rate = dot!(dtω, E_RT,f);
	store(heat_release_rate, rates, 0*stride, f);
	/*let Cc/*Cp|Cv*/ = a.iter().map(|a| dot(a[0], [(a[1], T), (a[2], T2), (a[3], T3), (a[4], T4)]));
	let m_rcp_ΣCCc = div(f.c(-1.), f.dot(concentrations.iter().copied().zip(Cc)), f);
	let dtT_T = mul(m_rcp_ΣCCc, heat_release_rate, f);
	store(mul(dtT_T, T, f), rates, 0*stride, f);*/
	//let R_S_Tdtn = mul(f.rcp(total_concentration), f.cdot(molar_mass[0..len-1].iter().map(|w| 1. - w/molar_mass[len-1]).zip(dtω.iter().copied()), None).unwrap(), f);
	//let dtS_S = f.add(R_S_Tdtn, dtT_T);
	//store(f.mul(dtS_S, variable), rates, 1*stride, f);
	let f = function.dfg;
	let mut w = String::new();
	for instruction in function.layout.block_insts(entry_block) {
		match f.inst_results(instruction) {
		 [] => (),
		 [result] => write!(w, "float {} = ", result)?,
		 _ => unimplemented!(),
		};
		use cranelift::codegen::ir::InstructionData::*;
		fn i32_from(offset: impl Into<i32>) -> i32 { offset.into() } // Workaround impl Into !From
		match &f[instruction] {
			UnaryIeee64{imm, ..} => write!(w, "{}", imm)?,
			Unary{opcode, arg} => write!(w, "{}({})", opcode, arg)?,
			Load{arg, offset, ..} => write!(w, "{}[{}]", arg, i32_from(*offset)/(std::mem::size_of::<f64>() as i32))?,
			Store{args, offset, ..} => write!(w, "{}[{}] = {}", args[1], i32_from(*offset)/(std::mem::size_of::<f64>() as i32), args[0])?,
			Call{func_ref, args, ..} => write!(w, "{:?}({})", Intrinsic::try_from(func_ref.as_u32()).unwrap(), args.first(&f.value_lists).unwrap())?,
			Binary{opcode, args} => write!(w, "{}({}, {})", opcode, args[0], args[1])?,
			_ => unimplemented!("{:?}", f[instruction])
		};
		write!(w, ";\n")?;
	}
	w
}
