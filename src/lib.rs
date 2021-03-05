#![feature(const_generics, const_evaluatable_checked, non_ascii_idents, type_ascription, once_cell, in_band_lifetimes, array_map, trait_alias, unboxed_closures, fn_traits)]
#![allow(incomplete_features, non_upper_case_globals, non_snake_case, confusable_idents, uncommon_codepoints)]
#![allow(unused_variables, dead_code)]
use num::log;

pub const K : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)

pub mod model;
use model::{Element, Troe};

#[derive(PartialEq, Debug, /*Eq*/)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure : f64 = 101325. / (K*NA); // 1 atm
	pub const T_split : f64 = 1000.;
	/*pub fn a(&self, T: f32) -> &[f32; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn specific_heat_capacity(&self, T: f32) -> f32 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn specific_enthalpy(&self, T: f32) -> f32 { let a = self.a(T); a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T } // /R
	pub fn specific_enthalpy_T(&self, T: f32) -> f32 { let a = self.a(T); a[5]/T+a[0]+a[1]/2.*T+a[2]/3.*T*T+a[3]/4.*T*T*T+a[4]/5.*T*T*T*T } // /RT
	pub fn specific_entropy(&self, T: f32) -> f32 { let a = self.a(T); a[6]+a[0]*log(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T } // /R*/
}

#[derive(Clone, Copy)] pub struct RateConstant {
	pub log_preexponential_factor: f64,
	pub temperature_exponent: f64,
	pub activation_temperature: f64
}

/*pub fn log_arrhenius(RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: RateConstant, T: f32) -> f32 {
	log_preexponential_factor + temperature_exponent*num::log(T) - activation_temperature*(1./T)
}*/

impl From<model::RateConstant> for RateConstant {
	fn from(model::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: model::RateConstant) -> Self {
		const J_per_cal: f64 = 4.184;
		Self{log_preexponential_factor: log(preexponential_factor)/*-temperature_exponent*log(K)*/, temperature_exponent, activation_temperature: activation_energy*J_per_cal/(K*NA)}
	}
}

use {std::{convert::TryInto, lazy::SyncLazy}, linear_map::LinearMap as Map};
static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, C: 12.011, N: 14.0067, O: 15.999, Ar: 39.95}").unwrap()
	.into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

#[derive(Debug)] pub struct Species {
	pub molar_mass: Box<[f64]>,
	pub thermodynamics: Box<[NASA7]>,
	diameter: Box<[f64]>,
	well_depth_J: Box<[f64]>,
	polarizability: Box<[f64]>,
	permanent_dipole_moment: Box<[f64]>,
	rotational_relaxation: Box<[f64]>,
	internal_degrees_of_freedom: Box<[f64]>,
	pub heat_capacity_ratio: Box<[f64]>,
}

pub enum ReactionModel {
	Elementary,
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

pub struct Model {
	pub species: Species,
	pub reactions: Box<[Reaction]>,
	//pub transport_polynomials: TransportPolynomials,
}

impl Model {
pub fn new(model::Model{species, reactions, ..}: model::Model) -> Self {
	let species: Box<[_]> = (species.into():Vec<_>).into();
	use std::ops::Deref;
	let species = species.deref();
	pub fn eval<T, U>(v: impl IntoIterator<Item=T>, f: impl Fn(T)->U) -> Box<[U]> { v.into_iter().map(f).collect() }
	/*let species = eval(species, |(k,specie):&(_,model::Specie)| {
		let mut specie = specie.clone();
		for T in specie.thermodynamic.temperature_ranges.iter_mut() { *T *= K; } // K->J
		for piece in specie.thermodynamic.pieces.iter_mut() { for (order, a) in piece[0..6].iter_mut().enumerate() { *a /= f64::powi(K, (1+order) as i32); } } // /K^n->/J^n
		(k, specie)
	});*/
	let species = species.deref();
	let molar_mass = eval(species, |(_,s)| s.composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
	let thermodynamics = eval(species, |(_, model::Specie{thermodynamic: model::NASA7{temperature_ranges, pieces},..})| match temperature_ranges[..] {
		[_,Tsplit,_] if Tsplit == NASA7::T_split => NASA7(pieces[..].try_into().unwrap()),
		[min, max] if min < NASA7::T_split && NASA7::T_split < max => NASA7([pieces[0]; 2]),
		ref ranges => panic!("{:?}", ranges),
	});
	let diameter = eval(species, |(_,s)| s.transport.diameter_Å*1e-10);
	let well_depth_J = eval(species, |(_,s)| s.transport.well_depth_K * K);
	use model::Geometry::*;
	let polarizability = eval(species, |(_,s)| if let Linear{polarizability_Å3,..}|Nonlinear{polarizability_Å3,..} = s.transport.geometry { polarizability_Å3*1e-30 } else { 0. });
	let permanent_dipole_moment = eval(species, |(_,s)|
		if let Nonlinear{permanent_dipole_moment_Debye,..} = s.transport.geometry { permanent_dipole_moment_Debye*Cm_per_Debye } else { 0. });
	let rotational_relaxation = eval(species, |(_,s)| if let Nonlinear{rotational_relaxation,..} = s.transport.geometry { rotational_relaxation } else { 0. });
	let internal_degrees_of_freedom = eval(species, |(_,s)| match s.transport.geometry { Atom => 0., Linear{..} => 1., Nonlinear{..} => 3./2. });
	let heat_capacity_ratio = eval(species, |(_,s)| {
		let f = match s.transport.geometry { Atom => 3., Linear{..} => 5., Nonlinear{..} => 6. };
		1. + 2./f
	});
	let species_names = eval(species, |(name,_)| *name);
	let species_names = species_names.deref();
	let species_composition = eval(species, |(_,s)| &s.composition);
	let species = Species{molar_mass, thermodynamics, diameter, well_depth_J, polarizability, permanent_dipole_moment, rotational_relaxation, internal_degrees_of_freedom, heat_capacity_ratio};
	//let transport_polynomials = species.transport_polynomials();
	let reactions = eval(reactions.into():Vec<_>, |model::Reaction{ref equation, rate_constant, model}| {
		for side in equation { for (specie, _) in side { assert!(species_names.contains(&specie), "{}", specie) } }
		let [reactants, products] = iter::vec::eval(equation, |e| eval(species_names, |&s| *e.get(s).unwrap_or(&0)));
		let net = eval(products.into_iter().zip(reactants.into_iter()).take(reactants.len()-1), |(&a, &b)| a as i8 - b as i8);
		{let mut net_composition = Map::new();
			for (s, &ν) in net.into_iter().enumerate() {
				for (element, &count) in species_composition[s] {
					if !net_composition.contains_key(&element) { net_composition.insert(element, 0); }
					*net_composition.get_mut(&element).unwrap() += ν as i8 * count as i8;
				}
			}
			for (_, &ν) in &net_composition { assert!(ν == 0, "{:?} {:?}", net_composition, equation); }
		}
		let [Σreactants, Σproducts] = [reactants.iter().sum(), products.iter().sum()];
		let Σnet = Σproducts as i8 - Σreactants as i8;
		let from = |efficiencies:Map<_,_>| eval(species_names, |&specie| *efficiencies.get(specie).unwrap_or(&1.));
		Reaction{
			reactants, products, net, Σreactants, Σproducts, Σnet,
			rate_constant: rate_constant.into(),
			model: {use model::ReactionModel::*; match model {
				Elementary => ReactionModel::Elementary,
				ThreeBody{efficiencies} => ReactionModel::ThreeBody{efficiencies: from(efficiencies)},
				PressureModification{efficiencies, k0} => ReactionModel::PressureModification{efficiencies: from(efficiencies), k0: k0.into()},
				Falloff{efficiencies, k0, troe} => ReactionModel::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe},
			}},
		}
	});
	Model{species, /*transport_polynomials,*/ reactions}
}
pub fn len(&self) -> usize { self.species.molar_mass.len() }
}

pub struct State {
    pub temperature: f64,
    pub pressure: f64,
    pub volume: f64,
    pub amounts: Box<[f64]>
}

pub struct Simulation<'t> {
	pub species_names: Box<[&'t str]>,
	pub state: State,
	pub time_step: f64,
}

impl Simulation<'t> {
	pub fn new(model::Model{species, state, time_step, ..}: &model::Model<'t>) -> ron::Result<Self> {
		let species_names = species.iter().map(|(name,_)| *name).collect():Box<_>;

		let model::State{temperature, pressure, volume, amount_proportions} = state;
		let pressure = pressure/NA;
		let temperature = *temperature;//*K; // K->J
		let amount = pressure * volume / (temperature * K);
		for (specie,_) in amount_proportions { assert!(species_names.contains(specie)); }
		let amount_proportions = species_names.iter().map(|specie| *amount_proportions.get(specie).unwrap_or(&0.)).collect():Box<_>;
		let amounts = amount_proportions.iter().map(|amount_proportion| amount/amount_proportions.iter().sum::<f64>() * amount_proportion).collect();

		Ok(Self{
			species_names,
			time_step: *time_step,
			state: State{temperature, pressure, volume: *volume, amounts}
		})
	}
}

#[cfg(test)] mod test;

#[derive(PartialEq, Eq)] pub enum Property { Pressure, Volume }
#[derive(Clone, Copy)] pub struct Constant<const CONSTANT: Property>(pub f32);
#[derive(derive_more::Deref, Debug)] pub struct StateVector<const CONSTANT: Property>(pub Box<[f32/*; T,P|V,[S-1]*/]>);
pub type Derivative<const CONSTANT: Property> = StateVector<CONSTANT>;

impl State {
	//pub fn constant<const CONSTANT: Property>(&Self{pressure, volume, ..}: &Self) -> f32 { // arbitrary_self_types
	pub fn constant<const CONSTANT: Property>(&self) -> Constant<CONSTANT> { let Self{pressure, volume, ..} = self;
		Constant(*{use Property::*; match CONSTANT {Pressure => pressure, Volume => volume}} as f32)
	}
}

use std::array::IntoIter;

impl<const CONSTANT: Property> From<&State> for StateVector<CONSTANT> {
	fn from(State{temperature, pressure, volume, amounts}: &State) -> Self {
		Self(IntoIter::new([*temperature as f32, *{use Property::*; match CONSTANT {Pressure => volume, Volume => pressure}} as f32]).chain(amounts[..amounts.len()-1].iter().map(|&n| n as f32)).collect())
	}
}

impl State {
	pub fn new<const CONSTANT: Property>(total_amount: f64, Constant(thermodynamic_state_constant): Constant<CONSTANT>, u: &StateVector<CONSTANT>) -> Self {
		let u = &u.0;
		let amounts = &u[2..];
		let (pressure, volume) = {use Property::*; match CONSTANT {
			Pressure => (thermodynamic_state_constant, u[1]),
			Volume => (u[1], thermodynamic_state_constant)
		}};
		State{
			temperature: u[0] as f64,
			pressure: pressure as f64,
			volume: volume as f64,
			amounts: amounts.iter().chain(&[total_amount as f32 - iter::into::Sum::<f32>::sum(amounts)]).map(|&n| n as f64).collect()
		}
	}
}

impl std::fmt::Display for State {
	fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
		let Self{temperature, pressure, volume, amounts} = self;
		write!(fmt, "T: {}, P: {}, V: {}, n: {:?}", temperature/*/K*/, pressure*NA, volume, amounts)
	}
}

use {cranelift::prelude::{*, types::{I32, F32}, codegen::{ir, binemit}}, cranelift_module::{Linkage, Module}};

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

struct Constants {
	_1: Value,
	_1_2: Value,
	_127: Value,
	exponent: Value,
	mantissa: Value,
	exp: [Value; 3],
	log: [Value; 3],
	constants: std::collections::HashMap<u32, Value>
}
impl Constants {
	fn new(f: &mut FunctionBuilder<'t>) -> Self { Self{
		_1: f![f f32const(1.)],
		_1_2: f![f f32const(1./2.)],
		_127: f![f iconst(I32, 127)],
		exponent: f![f iconst(I32, 0x7F800000)],
		mantissa: f![f iconst(I32, 0x007FFFFF)],
		exp: [1.0017247, 0.65763628, 0.33718944].map(|c| f![f f32const(c)]),
		log: [2.28330284476918490682, -1.04913055217340124191, 0.204446009836232697516].map(|c| f![f f32const(c)]),
		constants: Default::default(),
	}}
	#[track_caller] fn c(&mut self, f: &mut FunctionBuilder<'t>, constant: f32) -> Value {
		assert!(constant.is_finite(), "{}", constant);
		*self.constants.entry(constant.to_bits()).or_insert_with(|| f![f f32const(constant)])
	}
}

fn exp2(x: Value, Constants{_1_2, _127, exp, ..}: &Constants, f: &mut FunctionBuilder<'t>) -> Value {
	let ipart = f![f fcvt_to_sint(I32, f![f fsub(x, *_1_2)])];
	let fpart = f![f fsub(x, f![f fcvt_from_sint(F32, ipart)])];
	let expipart = f![f bitcast(F32, f![f ishl_imm(f![f iadd(ipart, *_127)], 23)])];
	let expfpart = fma![f (fma![f (exp[2], fpart, exp[1])], fpart, exp[0])];
	f![f fmul(expipart, expfpart)]
}

fn log2(x: Value, Constants{_1, exponent, _127, mantissa, log, ..}: &Constants, f: &mut FunctionBuilder<'t>) -> Value {
	let i = f![f /*raw_?*/bitcast(I32, x)];
	let e = f![f fcvt_from_sint(F32, f![f isub(f![f ushr_imm(f![f band(i, *exponent)], 23)], *_127)])];
	let m = f![f bor(f![f bitcast(F32, f![f band(i, *mantissa)])], *_1)];
	let p = fma![f (fma![f (log[2], m, log[1])], m, log[0])];
	let p = f![f fmul(p, f![f fsub(m, *_1)])]; //?
	f![f fadd(p, e)]
}

use std::f64::consts::LN_2;

fn log_arrhenius(RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: RateConstant,
														rcpT: Value, logT: Value, f: &mut FunctionBuilder<'t>) -> Value {
	let A = f![f f32const((log_preexponential_factor/LN_2) as f32)];
	let AeTβ = if temperature_exponent == 0. { A } else { fma![f (f![f f32const(temperature_exponent as f32)], logT, A)] };
	if activation_temperature == 0. { AeTβ } else { fma![f (f![f f32const((-activation_temperature/LN_2) as f32)], rcpT, AeTβ)] }
}

fn fdot<'t>(iter: impl IntoIterator<Item=(Value, impl FnOnce(&mut Constants, &mut FunctionBuilder<'t>)->Value)>, C: &mut Constants, f: &mut FunctionBuilder<'t>) -> Value {
	let mut iter = iter.into_iter();
	let mut sum = {let (a,b) = iter.next().unwrap(); f![f fmul(a, b(C, f))]};
	for (a,b) in iter { sum = fma![f (a, b(C, f), sum)]; }
	sum
}

#[track_caller] fn dot<T>(iter: impl IntoIterator<Item=(T, Value)>, mut sum: Option<Value>, C: &mut Constants, f: &mut FunctionBuilder<'t>) -> Option<Value>
where T: num::IsZero + num::IsOne + num::IsMinusOne + Into<f64> {
	for (c,v) in iter.into_iter() {
		//use num::{IsZero, IsOne, IsMinusOne};
		if c.is_zero() {}
		else if c.is_one() { sum = Some(match sum { Some(sum) => f![f fadd(sum, v)], None => v}); }
		else if c.is_minus_one() { sum = Some(match sum { Some(sum) => f![f fsub(sum, v)], None => v}); }
		else { let c = C.c(f, c.into() as f32); sum = Some(match sum { Some(sum) => fma![f (c, v, sum)], None => f![f fmul(c, v)] }); }
	}
	sum
}


impl ReactionModel {
fn efficiency(&self, f: &mut FunctionBuilder<'t>, C: &mut Constants, logT: Value, rcpT: Value, mT: Value, mrcpT: Value, concentrations: &[Value], log_k_inf: Value) -> Value {
	use ReactionModel::*; match self {
		Elementary => C._1,
		ThreeBody{efficiencies} => { dot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None, C, f).unwrap() },
		PressureModification{efficiencies, k0} => {
			let Pr = f![f fmul(dot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None, C, f).unwrap(), exp2(f![f fsub(log_arrhenius(*k0, rcpT, logT, f), log_k_inf)], C, f))];
			f![f fdiv(Pr, f![f fadd(C._1, Pr)])]
		}
		Falloff{efficiencies, k0, troe} => {
			let Pr = f![f fmul(dot(efficiencies.iter().copied().zip(concentrations.iter().copied()), None, C, f).unwrap(), exp2(f![f fsub(log_arrhenius(*k0, rcpT, logT, f), log_k_inf)], C, f))];
			let model::Troe{A, T3, T1, T2} = *troe;
			let Fcent = fma![f (f![f f32const(1.-A as f32)], exp2(f![f fdiv(mT, f![f f32const(T3 as f32)])], C, f),
														 fma![f (f![f f32const(A as f32)], exp2(f![f fdiv(mT, f![f f32const(T1 as f32)])], C, f),
																																							exp2(f![f fmul(mrcpT, f![f f32const(T2 as f32)])], C, f) )] )];
			let logFcent = log2(Fcent, C, f);
			let c =fma![f (C.c(f, -0.67), logFcent, C.c(f, -0.4*f64::log2(10.) as f32))];
			let N = fma![f (C.c(f, -1.27), logFcent, C.c(f, 0.75*f64::log2(10.) as f32))];
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

impl Model {
pub fn rate<const CONSTANT: Property>(&self) -> (Box<[Trap]>, (extern fn(f32, *const f32, *mut f32), /*usize*/()), impl Fn(Constant<CONSTANT>, &StateVector<CONSTANT>, &mut Derivative<CONSTANT>)/*+'static*/) {
	let builder = cranelift_jit::JITBuilder::new(cranelift_module::default_libcall_names());
	let mut module = cranelift_jit::JITModule::new(builder);
  let mut context = module.make_context();
	let PTR = module.target_config().pointer_type();
  let params = [("constant", F32), ("state", PTR), ("rate", PTR)];
	context.func.signature.params = params.iter().map(|(_,r#type)| AbiParam::new(*r#type)).collect();
	let mut function_builder_context = FunctionBuilderContext::new();
	let mut builder = FunctionBuilder::new(&mut context.func, &mut function_builder_context);
	let entry_block = builder.create_block();
	builder.append_block_params_for_function_params(entry_block);
	builder.switch_to_block(entry_block);
	builder.seal_block(entry_block);
	let [constant, state, rate]: [Value; 3] = builder.block_params(entry_block).try_into().unwrap();
	let flags = MemFlags::new();
	let ref mut f = builder;
	let ref mut C = Constants::new(f);
	let _m1 = f![f f32const(-1.)];
	use std::mem::size_of;
	let T = f![f load(F32, flags, state, 0*size_of::<f32>() as i32)];
	let logT = log2(T, C, f);
	let rcpT = f![f fdiv(C._1, T)];
	let T2 = f![f fmul(T, T)];
	let T3 = f![f fmul(T2, T)];
	let T4 = f![f fmul(T3, T)];
	let mT = f![f fneg(T)];
	let mrcpT = f![f fneg(rcpT)];
	let Self{species: Species{molar_mass: W, thermodynamics, heat_capacity_ratio, ..}, reactions} = self;
	let len = self.len();
	let a = thermodynamics.iter().map(|s| s.0[1]).collect(): Box<_>;
	let G_RT =
		a.iter().map(|a| dot(
			IntoIter::new([(a[5], rcpT), (a[0]*LN_2, logT), (-a[1]/2., T), ((1./3.-1./2.)*a[2], T2), ((1./4.-1./3.)*a[3], T3), ((1./5.-1./4.)*a[4]/5., T4)]),
			Some(f![f f32const((a[0]-a[6]) as f32)]), C, f).unwrap()).collect(): Box<_>;
	let logP0_RT = f![f fsub(f![f f32const(f64::log2(NASA7::reference_pressure) as f32)], logT)];
	let variable = f![f load(F32, flags, state, 1*size_of::<f32>() as i32)];
	let (pressure, volume) = {use Property::*; match CONSTANT {Pressure => (constant, variable), Volume => (variable, constant)}};
	let kT = f![f fmul(C.c(f, K as f32), T)];
	let total_concentration = f![f fdiv(pressure/*/Na*/, kT)]; // n/V = P/RT
	let amounts = (0..len-1).map(|i| f![f load(F32, flags, state, ((2+i)*size_of::<f32>()) as i32)]).collect(): Box<_>;
	let rcpV = f![f fdiv(C._1, volume)];
	let concentrations = amounts.iter().map(|&n| f![f fmul(n/*.max(0.)*/, rcpV)]).collect(): Box<_>;
	let Ca = f![f fsub(total_concentration, dot(std::iter::repeat(1.).zip(concentrations.iter().copied()), None, C, f).unwrap())];
	//if Ca < 0. { dbg!(T, C, concentrations, Ca); throw!(); }
	let ref concentrations = [&concentrations as &[_],&[Ca]].concat();
	let log_concentrations = concentrations.iter().map(|&x| log2(x, C, f)).collect(): Box<[Value]>;
	let mut dtω = (0..len-1).map(|_| None).collect(): Box<_>;
	for Reaction{reactants, products, net, Σnet, rate_constant, model, ..} in reactions.iter() {
		let log_k_inf = log_arrhenius(*rate_constant, rcpT, logT, f);
		let c = model.efficiency(f, C, logT, rcpT, mT, mrcpT, concentrations, log_k_inf); // todo: CSE
		let Rf = exp2(dot(reactants.iter().copied().zip(log_concentrations.iter().copied()), Some(log_k_inf), C, f).unwrap(), C, f);
		let m_log_equilibrium_constant = dot(net.iter().copied().chain(IntoIter::new([-Σnet])).zip(G_RT.iter().chain(&[logP0_RT]).copied()), None, C, f).unwrap();
		let Rr = exp2(dot(products.iter().copied().zip(log_concentrations.iter().copied()), Some(f![f fadd(log_k_inf, m_log_equilibrium_constant)]), C, f).unwrap(), C, f);
		let R = f![f fsub(Rf, Rr)];
		let cR = f![f fmul(c, R)];
		for (index, &ν) in net.iter().enumerate() {
			let dtω = &mut dtω[index];
			match ν {
					0 => {},
					1 => match dtω { None => *dtω = Some(cR), Some(dtω) => *dtω = f![f fadd(*dtω, cR)] }
					-1 => match dtω { None => *dtω = Some(f![f fneg(cR)]), Some(dtω) => *dtω = f![f fsub(*dtω, cR)] }
					ν => match dtω { None => *dtω = Some(f![f fmul(C.c(f, ν as f32), cR)]), Some(dtω) => *dtω = fma![f (C.c(f, ν as f32), cR, *dtω)] }
			}
		}
	}
	//for (specie, dtω) in dtω.iter().enumerate() { if let None = dtω { panic!("Inactive specie {}", specie); } }

	let a = {use Property::*; match CONSTANT {Pressure => a, Volume => a.iter().zip(heat_capacity_ratio.iter()).map(|(a,γ)| a.map(|a| a / γ)).collect()}};
	struct Dot<const N: usize>(f64, [(f64, Value); N]);
	impl<'t, const N: usize> FnOnce<(&mut Constants, &mut FunctionBuilder<'t>,)> for Dot<N> {
		type Output = Value;
		extern "rust-call" fn call_once(self, (C, f,): (&mut Constants, &mut FunctionBuilder<'t>,)) -> Self::Output {
			dot(IntoIter::new(self.1), Some(f![f f32const(self.0 as f32)]), C, f).unwrap()
		}
	}
	let Cc/*Cp|Cv*/ = a.iter().map(|a| Dot(a[0], [(a[1], T), (a[2], T2), (a[3], T3), (a[4], T4)]));
	let m_rcp_ΣCCc = f![f fdiv(_m1, fdot(concentrations.iter().copied().zip(Cc), C, f))];
	let E_T/*H/T|U/T*/ = a.iter().map(|a| Dot(a[0], [(a[5], rcpT), (a[1]/2., T), (a[2]/3., T2), (a[3]/4., T3), (a[4]/5., T4)]));
	let dtω = dtω.into_iter().map(|dtω| dtω.unwrap()).collect(): Box<_>;
	let dtT_T = f![f fmul(m_rcp_ΣCCc, fdot(dtω.iter().copied().zip(E_T), C, f))];
	f![f store(flags, f![f fmul(dtT_T, T)], rate, 0*size_of::<f32>() as i32)];
	let R_S_Tdtn = f![f fmul(f![f fdiv(kT, pressure/*/Na*/)], dot(W[0..len-1].iter().map(|w| 1. - w/W[len-1]).zip(dtω.iter().copied()), None, C, f).unwrap())]; // R/A Tdtn (constant pressure: A=V, constant volume: A=P)
	let dtS_S = f![f fadd(R_S_Tdtn, dtT_T)];
	f![f store(flags, f![f fmul(dtS_S, variable)], rate, 1*size_of::<f32>() as i32)];
	let dtn = dtω.into_iter().map(|&dtω| f![f fmul(volume, dtω)]).collect(): Box<_>;
	for (i, &dtn) in dtn.into_iter().enumerate() { f![f store(flags, dtn, rate, ((2+i)*size_of::<f32>()) as i32)]; }
	builder.ins().return_(&[]);
	builder.finalize();
	let clif = builder.display(None);
	//eprintln!("{}", clif);
	/*{
		let function = cranelift_reader::parse_functions(&clif.to_string()).unwrap().remove(0);
		let mut context = codegen::Context::new();
		context.func = function;
		let mut mem = vec![];
		let isa_builder = isa::lookup(target_lexicon::Triple::host()).unwrap();
		let mut flag_builder = settings::builder();
		flag_builder.enable("is_pic").unwrap();
		let isa = isa_builder.finish(settings::Flags::new(flag_builder));
		let code_info = context.compile_and_emit(&*isa, &mut mem, &mut binemit::NullRelocSink{}, &mut binemit::NullTrapSink{}, &mut binemit::NullStackMapSink{});
		use capstone::arch::BuildsCapstone;
		let capstone = capstone::Capstone::new().x86().mode(capstone::arch::x86::ArchMode::Mode64).build().unwrap();
		let instructions = capstone.disasm_all(&mem, 0).unwrap();
		//for i in instructions.iter() { println!("{}\t{}", i.mnemonic().unwrap(), i.op_str().unwrap()); }
	}*/
  let id = module.declare_function(&"", Linkage::Export, &context.func.signature).unwrap();
  struct Traps (Vec<Trap>);
  impl Traps { fn new() -> Self { Self(Vec::new()) } }
	impl binemit::TrapSink for Traps {
    fn trap(&mut self, code_offset: CodeOffset, source_location: SourceLoc, trap_code: TrapCode) {
			self.0.push(Trap{code_offset, source_location, trap_code});
    }
	}
	let mut trap_sink = Traps::new();
	module.define_function(id, &mut context, &mut trap_sink).unwrap();
	module.finalize_definitions();
	//let (function, size) = module.get_finalized_function(id);
	let (function, size) = (module.get_finalized_function(id), ());
	let function = unsafe{std::mem::transmute::<_,extern fn(f32, *const f32, *mut f32)>(function)};
	(trap_sink.0.into_boxed_slice(), (function, size), move |constant, state, derivative| {
		let constant = constant.0 as f32;
		function(constant, state.0.as_ptr(), derivative.0.as_mut_ptr());
	})
}
}
