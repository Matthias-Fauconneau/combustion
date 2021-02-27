#![feature(const_generics, const_evaluatable_checked, non_ascii_idents, type_ascription, once_cell)]
#![allow(incomplete_features, non_upper_case_globals, non_snake_case, confusable_idents)]
use std::ops::Deref;
use num::log;

pub const K : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)

use model::{Element, Troe};

#[derive(PartialEq, Debug, /*Eq*/)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure : f64 = 101325. / NA; // 1 atm
	/*pub const T_split : f64 = 1000.*K;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn specific_enthalpy(&self, T: f64) -> f64 { let a = self.a(T); a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T } // /R
	pub fn specific_enthalpy_T(&self, T: f64) -> f64 { let a = self.a(T); a[5]/T+a[0]+a[1]/2.*T+a[2]/3.*T*T+a[3]/4.*T*T*T+a[4]/5.*T*T*T*T } // /RT
	pub fn specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); a[6]+a[0]*log(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T } // /R*/
}

#[derive(Clone, Copy)] pub struct RateConstant {
	pub log_preexponential_factor: f64,
	pub temperature_exponent: f64,
	pub activation_temperature: f64
}

/*pub fn log_arrhenius(RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: RateConstant, T: f64) -> f64 {
	log_preexponential_factor + temperature_exponent*num::log(T) - activation_temperature*(1./T)
}*/

impl From<model::RateConstant> for RateConstant {
	fn from(model::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: model::RateConstant) -> Self {
		const J_per_cal: f64 = 4.184;
		Self{log_preexponential_factor: log(preexponential_factor)-temperature_exponent*log(K), temperature_exponent, activation_temperature: activation_energy*J_per_cal/NA}
	}
}

#[derive(Clone, Copy, Debug, PartialEq)] pub struct State<const S: usize> {
    pub temperature: f64,
    pub pressure: f64,
    pub volume: f64,
    pub amounts: [f64; S]
}

use std::{convert::TryInto, lazy::SyncLazy, linear_map::LinearMap as Map};
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

#[derive(Debug)] pub enum ReactionModel {
	Elementary,
	ThreeBody { efficiencies: Box<[f64]> },
	PressureModification { efficiencies: Box<[f64]>, k0: RateConstant },
	Falloff { efficiencies: Box<[f64]>, k0: RateConstant, troe: Troe },
}

#[derive(Debug)] pub struct Reaction {
	pub reactants: Box<[u8]>,
	pub products: Box<[u8]>,
	pub net: Box<[i8/*; S-1*/]>,
	pub Σreactants: u8,
	pub Σproducts: u8,
	pub Σnet: i8,
	pub rate_constant: RateConstant,
	pub model: ReactionModel,
}

#[derive(Debug)] pub struct Model {
	pub species: Species,
	pub reactions: Box<[Reaction]>,
	//pub transport_polynomials: TransportPolynomials,
}

impl Model {
pub fn new(model::Model{species, reactions, ..}: model::Model) -> Self {
	let species: Box<[_]> = (species.into():Vec<_>).into();
	let species = species.deref();
	pub fn eval<T, U>(v: impl IntoIterator<Item=T>, f: impl Fn(T)->U) -> Box<[U]> { v.into_iter().map(f).collect() }
	let species = eval(species, |(k,specie):&(_,model::Specie)| {
		let mut specie = specie.clone();
		for T in specie.thermodynamic.temperature_ranges.iter_mut() { *T *= K; } // K->J
		for piece in specie.thermodynamic.pieces.iter_mut() { for (order, a) in piece[0..6].iter_mut().enumerate() { *a /= f64::powi(K, (1+order) as i32); } } // /K^n->/J^n
		(k, specie)
	});
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
	let species = Species_{molar_mass, thermodynamics, diameter, well_depth_J, polarizability, permanent_dipole_moment, rotational_relaxation, internal_degrees_of_freedom, heat_capacity_ratio};
	//let transport_polynomials = species.transport_polynomials();
	let reactions = eval(reactions.into():Vec<_>, |model::ReactionModel{ref equation, rate_constant, model}| {
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
				Falloff{efficiencies, k0, troe: Troe{A,T3,T1,T2}} => ReactionModel::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe: Troe{A,T3:T3*K,T1:T1*K,T2:T2*K}},
			}},
		}
	});
	System{species, /*transport_polynomials,*/ reactions}
}
pub fn len(&self) -> usize { self.species.molar_mass.len() }
}

pub struct Simulation<'t, const S: usize> where [(); S-1]: {
	pub species_names: [&'t str; S],
	pub state: State<S>,
	pub time_step: f64,
}

impl<const S: usize /*= SPECIES_LEN*/> Simulation<'t, S> where [(); S-1]: {
	pub fn new(simulation: &'t str) -> ron::Result<Self> where /*'b: 't,*/ [(); S]: {
		let model::Model{species, state, time_step, ..} = ::ron::de::from_str(simulation)?;
		let species: Box<[_]> = (species.into():Vec<_>).into();
		let species : Box<[_; S]> = {let len = species.len(); unwrap::unwrap!(species.try_into(), "Compiled for {} species, got {}", S, len)};
		let species = species.deref();
		let species_names = eval(species, |(name,_)| *name);

		let model::State{temperature, pressure, volume, amount_proportions} = state;
		let pressure = pressure/NA;
		let temperature = temperature*K; // K->J
		let amount = pressure * volume / temperature;
		for (specie,_) in &amount_proportions { assert!(species_names.contains(specie)); }
		let amount_proportions = eval(species_names, |specie| *amount_proportions.get(specie).unwrap_or(&0.));
		let amounts = eval(amount_proportions/*.prefix()*/, |amount_proportion| amount/amount_proportions.iter().sum::<f64>() * amount_proportion);

		Ok(Self{
			species_names,
			time_step,
			state: State{temperature, pressure, volume, amounts}
		})
	}
}

#[cfg(test)] mod test;

#[derive(PartialEq, Eq)] pub enum Property { Pressure, Volume }
pub struct Constant<const CONSTANT: Property>(f64);
#[derive(Debug, derive_more::Deref)] pub struct StateVector<const CONSTANT: Property>(pub Box<[f64/*; T,P|V,[S-1]*/]>);
pub type Derivative<const CONSTANT: Property> = StateVector<CONSTANT>;

impl State {
	//pub fn constant<const CONSTANT: Property>(&Self{pressure, volume, ..}: &Self) -> f64 { // arbitrary_self_types
	pub fn constant<const CONSTANT: Property>(&self) -> Constant<CONSTANT> { let Self{pressure, volume, ..} = self;
		Constant(*{use Property::*; match CONSTANT {Pressure => pressure, Volume => volume}})
	}
}

impl<const CONSTANT: Property> From<&State> for StateVector<CONSTANT> {
	fn from(State{temperature, pressure, volume, amounts}: &State) -> Self {
		Self([*temperature, *{use Property::*; match CONSTANT {Pressure => volume, Volume => pressure}}], amounts[..amounts.len()-1]])
	}
}

impl State {
	pub fn new<const CONSTANT: Property>(total_amount: f64, Constant(thermodynamic_state_constant): Constant<CONSTANT>, u: &State<CONSTANT, S>) -> Self {
		let u = u.0;
		let amounts: &[_; S-1] = u.suffix();
		let (pressure, volume) = {use Property::*; match CONSTANT {
			Pressure => (thermodynamic_state_constant, u[1]),
			Volume => (u[1], thermodynamic_state_constant)
		}};
		State{temperature: u[0], pressure, volume, amounts: from_iter(amounts.copied().chain([total_amount - iter::into::Sum::<f64>::sum(amounts)]))}
	}
}

impl<const S: usize> std::fmt::Display for State<S> {
	fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
		let Self{temperature, pressure, volume, amounts} = self;
		write!(fmt, "T: {}, P: {}, V: {}, n: {:?}", temperature/K, pressure*NA, volume, amounts)
	}
}

use crate::frontend::*;
use cranelift::prelude::*;
use cranelift_jit::{JITBuilder, JITModule};
use cranelift_module::{DataContext, Linkage, Module};

/*#![feature(non_ascii_idents, type_ascription, once_cell, array_map, proc_macro_quote, in_band_lifetimes)]
#![allow(confusable_idents, non_upper_case_globals, non_snake_case, unused_variables)]
use std::ops::Deref;
use system::{Model, Reaction, System};
use proc_macro::TokenStream;
fn literal(value: impl std::fmt::Debug) -> TokenStream { format!("{:?}", value).parse().unwrap() }
//use syn::{ItemFn, Expr, parse_quote as quote};
use proc_macro::quote; type Expr = TokenStream; type ItemFn = TokenStream;*/

fn exp2(f: FunctionBuilder, x: Value) {
	use cranelift_codegen::ir::types::{F32, F64};
	let f = || f.ins();
	let x = f().fdemote(F32, x);
	let ipart = f().fcvt_to_sint(f().fsub(x, f().f32const(1./2.)));
	let fpart = f().sub(x, f().fcvt_from_sint(ipart));
	let expipart = f().bitcast(f().ushl_imm(f().add(ipart, f.iconst(127)), 23));
	let c = [1.0017247, 0.65763628, 0.33718944].map(|c| f().f32const(c));
	let expfpart = f().fma(f().fma(c[2], x, c[1]), x, c[0]);
	f().fpromote(F64, f().fmul(expipart, expfpart))
}

fn log2(f: FunctionBuilder, x: Value) {
	use cranelift_codegen::ir::types::{F32, F64};
	let f = || f.ins();
	let x = f().fdemote(F32, x);
	let exponent = f().iconst(0x7F800000);
	let mantissa = f().iconst(0x007FFFFF);
	let i = f()./*raw_?*/bitcast(x);
	let e = f().fcvt_from_sint(f().sub(ins().ushr_imm(f().and(i, exponent), 23), f().iconst(127)));
	let m = f().or(f().bitcast(f().and(i, mantissa)), _1);
	let c = [2.28330284476918490682, -1.04913055217340124191, 0.204446009836232697516].map(|c| f().f32const(c));
	let p = f().fma(f().fma(c[2], x, c[1]), x, c[0]);
	let p = f().fmul(p, f().fsub(m, _1)); //?
	f().fpromote(F64, f().fadd(p, e))
}

use std::f32::consts::LN_2;

fn log_arrhenius(RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: RateConstant, f: FunctionBuilder,
														rcpT: Value, logT: Value) -> Value {
	let f = || f.ins();
	f().fma(f().f64const(-activation_temperature/LN_2), rcpT, f().fma(f().f64const(temperature_exponent), logT, f().f64const(log_preexponential_factor/LN_2)))
}

fn dot<T: num::IsZero + num::IsOne + num::IsMinusOne>(f: FunctionBuilder, zip: impl IntoIterator<Item=(T, Value)) -> Value {
	let mut sum = {let (a,b) = zip.next().unwrap(); f().fmul(a, b)};
	for (a,b) in zip { sum = f().fma(c, v, sum); }
	sum
}

fn fdot(f: FunctionBuilder, zip: impl IntoIterator<Item=(Value, Value), mut sum: Option<Value>) -> Option<Value> {
	for (c,v) in zip {
		if c.is_zero() {}
		else if c.is_one() { sum = Some(match sum { Some(sum) => f().fadd(sum, v), None => v}); }
		else if c.is_minus_one() { sum = Some(match sum { Some(sum) => f().fsub(sum, v), None => v}); }
		else { let c = f().f32const(c); sum = Some(match sum { Some(sum) => f().fma(c, v, sum), None => f().fmul(c, v) }); }
	}
	sum
}


impl Model {
fn efficiency(&self, &mut f: FunctionBuilder, _1: Value, logT: Value, rcpT: Value, mT: Value, concentrations: &[Value], log_k_inf: Value) -> Value {
	use Model::*;
	match model {
		Elementary => _1,
		ThreeBody{efficiencies} => { dot(f, efficiencies, concentrations, None) },
		PressureModification{efficiencies, k0} => {
			let Pr = f().fmul(dot(efficiencies, concentrations, None), exp2(f, f().fsub(log_arrhenius(f, k0, rcpT, logT), log_k_inf))); // [k0/kinf] = [1/C] (m3/mol)
			let f().fdiv(Pr, f().add(_1, Pr))
		}
		Falloff{efficiencies, k0, troe} => {
			let Pr = f().mul(dot(efficiencies, concentrations, None), exp(f, f().sub(log_arrhenius(f, k0, rcpT, logT), log_k_inf))); // [k0/kinf] = [1/C] (m3/mol)
			let model::Troe{A, T3, T1, T2} = troe;
			let logFcent = log2(f, f().fma(f().f64const(1.-A), exp(f, f().fdiv(mT/T3)), f.fma(f().f64const(A), exp(f,f().fdiv(mT, T1)), exp(f, f().fmul(mrcpT, T2)))));
			let C =f().fma(f().f64const(-0.67), logFcent, f().f64const(-0.4*f64::log2(10.)));
			let N = f().fma(f().f64const(-1.27), logFcent, f().f64const(0.75*f64::log2(10.)));
			let logPr᛭C = f().add(log2(f, Pr), C);
			let f1 = f().fdiv(logPr᛭C, f().fma(-0.14, logPr᛭C, N));
			let F = exp2(f, f().fdiv(logFcent, f().fma(f1, f1, _1)));
			f().fdiv(Pr, f().fma(Pr, F, F))
		}
	}
}
}

impl Model {
pub fn rate<const CONSTANT: Property>(&self) -> fn(Constant<CONSTANT>, State<CONSTANT>, &mut Derivative<CONSTANT>) {
	let builder = JITBuilder::new(cranelift_module::default_libcall_names());
	let builder_context: FunctionBuilderContext::new(),
	let contex =: module.make_context();
	let data_context = DataContext::new(),
  let module = JITModule::new(builder);
  use cranelift_codegen::ir::types::F64;
  let PTR= module.target_config().pointer_type();
  let params = [("constant", F64), ("state", PTR), ("rate", PTR)];
	context.func.signature.params = params.map(|(_,type)| AbiParam::new(type)).collect();
	let mut builder = FunctionBuilder::new(&mut context.func, &mut builder_context);
	let entry_block = builder.create_block();
	builder.append_block_params_for_function_params(entry_block);
	builder.switch_to_block(entry_block);
	builder.seal_block(entry_block);
	let mut variables = HashMap::new();
	let declare_variable = {
		let mut index = 0;
		|builder: &mut FunctionBuilder, /*variables: &mut HashMap<String, Variable>,*/ name: &str, type: Type| -> Variable {
			let variable = Variable::new(*index);
			if !variables.contains_key(name) {
					variables.insert(name.into(), variable);
					builder.declare_var(variable, type);
					index += 1;
			}
			variable
		}
	};
	let [constant, state, rate] = params.iter().enumerate().map(|(i, (name, type))|{
		let var = declare_variable(builder, name, type); builder.def_var(var, builder.block_params(entry_block)[i]); var
	}).collect();
	let flags = MemFlags::new();
	let f = || builder.ins();
	let _1 = ins().f64const(1.);
	let _m1 = ins().f64const(-1.);
	use std::mem::size_of;
	let T = f().load(F64, flags, state, 0*size_of::<f64>())
	let logT = log(ins, T);
	let rcpT = f().fdiv(_1, T);
	let T2 = f().fmul(T, T);
	let T3 = f().fmul(T2, T);
	let T4 = f().fmul(T3, T);
	let mT =f().fneg(T);
	let mrcpT = f().fneg(mT);
	let Self{species: Species{molar_mass: W, thermodynamics}, reactions}
	let a = thermodynamics.iter().map(|s| s.0[1]).collect();
	let G_RT =
		a.iter().map(|a|
			dot(f, [(a[5], Tr), (a[0]*LN_2, logT), (-a[1]/2., T), ((1./3.-1./2)*a[2], T2), ((1./4.-1./3)*a[3], T3), ((1./5.-1./4)*a[4]/5., T4)], f().f32const(a[0]-a[6]))
		).collect();
	let variable = f().load(F64, flags, state, 1*size_of::<f64>());
	let (pressure, volume) = {use Property::*; match CONSTANT {Pressure => (constant, variable), Volume => (variable, constant)}};
	let C = pressure / T; // n/V = P/kT
	let len = model.len();
	let amounts = (0..len-1).map(|i| f().load(F64, flags, state, (2+i)*size_of::<f64>())).collect();
	let rcpV = f().fdiv(_1, volume);
	let concentrations = amounts.iter().map(|n| f().fmul(n/*.max(0.)*/, rcpV)).collect();
	let Ca = f().fsub(C, dot(f, std::iter::repeat(1.).zip(concentrations), None));
	//if Ca < 0. { dbg!(T, C, concentrations, Ca); throw!(); }
	let ref concentrations = [concentrations,[Ca]].concat();
	let ref log_concentrations = concentrations.iter().map(|x| log2(f, x)).collect();
	let mut dtω = (0..len-1).map(|_| None).collect();
	for Reaction{reactants, products, net, Σnet, rate_constant, model, ..} in self.reactions {
		let log_k_inf = log_arrhenius(rate_constant, rcpT, logT);
		let c = model.efficiency(logT, rcpT, mT, concentrations, log_k_inf); // todo: CSE
		let Rf = exp2(f, dot(reactants.iter().zip(log_concentrations), Some(log_k_inf)));
		let m_log_equilibrium_constant = dot(net.iter().zip(G_RT), Some(Σnet*logP0_RT));
		let Rr = exp2(f, dot(products.iter().zip(log_concentrations), Some(log_k_inf + m_log_equilibrium_constant)));
		//assert!(Rr.is_finite(), "{} {} {}", $dot_products(log_concentrations), log_k_inf, log_equilibrium_constant);
		let R = f().fsub(Rf, Rr);
		let cR = f().fmul(c * R);
		for (index, ν) in net.iter().enumerate() {
			let ref mut dtω = dtω[index];
			match ν {
					0 => {},
					1 => match dtω { None => dtω = Some(cR), Some(dtω) => dtω = f().fadd(dtω, cR) }
					-1 => match dtω { None => dtω = Some(cR), Some(dtω) => dtω = f().fsub(dtω, cR) }
					ν => match dtω { None => dtω = Some(f().fmul(ν, cR)), Some(dtω) => dtω = f().fma(ν, cR, dtω) }
			}
		}
	}
	let Cc/*Cp|Cv*/ = a.iter().zip(heat_capacity_ratio).map(|(a,γ)|
		let a = {use Property::*; match CONSTANT {Pressure => a, Volume => eval(a, |a| a / γ)}};
		dot(f, [(a[1], T), (a[2], T2), (a[3], T3), (a[4]., T4)], f().f32const(a[0]));
	).collect();
	let m_rcp_ΣCCc = f().fdiv(_m1, fdot(f, concentrations.iter().zip(Cc)));
	let E_T/*H/T|U/T*/ = a.iter().zip(heat_capacity_ratio).map(|(a,γ)|
		let a = {use Property::*; match CONSTANT {Pressure => a, Volume => eval(a, |a| a / γ)}};
		dot(f, [(a[5], rcpT), (a[1]/2., T), (a[2]/3., T2), (a[3]/4., T3), (a[4]/5., T4)], f().f32const(a[0]));
	).collect();
	let dtT_T = m_rcp_ΣCCc * fdot(f, dtω.iter().zip(E_T));
	f().store(F64, flags, state, f().fmul(dtT_T, T), 0*size_of::<f64>());
	let R_S_Tdtn = f().fdiv(T, pressure) * dot(W[0..len-1].iter().map(|w| 1. - w/W[len-1]).zip(dtω), None); // R/A Tdtn (constant pressure: A=V, constant volume: A=P)
	let dtS_S = f().add(R_S_Tdtn, dtT_T);
	let dtn = dtω.map(|dtω| volume * dtω).collect();
	for (i, dtn) in dtn { f().store(F64, flags, state, n, (2+i)*size_of::<f64>());
	builder.finalize();
  let id = module.declare_function(&"", Linkage::Export, &context.func.signature)?;
	module.define_function(id, &mut context, &mut codegen::binemit::NullTrapSink{})?;
	module.finalize_definitions();
	module.get_finalized_function(id) as  fn(Constant<CONSTANT>, State<CONSTANT>) -> Derivative<CONSTANT>;
}
}
