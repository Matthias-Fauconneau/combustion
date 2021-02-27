#![feature(const_generics, const_evaluatable_checked, non_ascii_idents, type_ascription, once_cell, in_band_lifetimes, array_map, trait_alias)]
#![allow(incomplete_features, non_upper_case_globals, non_snake_case, confusable_idents, uncommon_codepoints)]
use std::ops::Deref;
use num::log;

pub const K : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)

mod model;
use model::{Element, Troe};

#[derive(PartialEq, Debug, /*Eq*/)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure : f64 = 101325. / NA; // 1 atm
	pub const T_split : f64 = 1000.*K;
	/*pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
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
				Falloff{efficiencies, k0, troe: Troe{A,T3,T1,T2}} => ReactionModel::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe: Troe{A,T3:T3*K,T1:T1*K,T2:T2*K}},
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
	pub fn new(simulation: &'t str) -> ron::Result<Self> {
		let model::Model{species, state, time_step, ..} = ::ron::de::from_str(simulation)?;
		let species: Box<[_]> = (species.into():Vec<_>).into();
		let species_names = species.iter().map(|(name,_)| *name).collect():Box<_>;

		let model::State{temperature, pressure, volume, amount_proportions} = state;
		let pressure = pressure/NA;
		let temperature = temperature*K; // K->J
		let amount = pressure * volume / temperature;
		for (specie,_) in &amount_proportions { assert!(species_names.contains(specie)); }
		let amount_proportions = species_names.iter().map(|specie| *amount_proportions.get(specie).unwrap_or(&0.)).collect():Box<_>;
		let amounts = amount_proportions.iter().map(|amount_proportion| amount/amount_proportions.iter().sum::<f64>() * amount_proportion).collect();

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
		Self([&[*temperature, *{use Property::*; match CONSTANT {Pressure => volume, Volume => pressure}}], &amounts[..amounts.len()-1]].concat().into_boxed_slice())
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
		State{temperature: u[0], pressure, volume, amounts: [amounts,&[total_amount - iter::into::Sum::<f64>::sum(amounts)]].concat().into_boxed_slice()}
	}
}

impl std::fmt::Display for State {
	fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
		let Self{temperature, pressure, volume, amounts} = self;
		write!(fmt, "T: {}, P: {}, V: {}, n: {:?}", temperature/K, pressure*NA, volume, amounts)
	}
}

use cranelift::prelude::{*, types::{I32, F32, F64}};
use cranelift_jit::{JITBuilder, JITModule};
use cranelift_module::{Linkage, Module};

fn exp2(f: &mut FunctionBuilder<'t>, x: Value) -> Value {
	let x = f.ins().fdemote(F32, x);
	let ipart = f.ins().fcvt_to_sint(I32, f.ins().fsub(x, f.ins().f32const(1./2.)));
	let fpart = f.ins().fsub(x, f.ins().fcvt_from_sint(F32, ipart));
	let expipart = f.ins().bitcast(F32, f.ins().ishl_imm(f.ins().iadd(ipart, f.ins().iconst(I32, 127)), 23));
	let c = [1.0017247, 0.65763628, 0.33718944].map(|c| f.ins().f32const(c));
	let expfpart = f.ins().fma(f.ins().fma(c[2], fpart, c[1]), fpart, c[0]);
	f.ins().fpromote(F64, f.ins().fmul(expipart, expfpart))
}

fn log2(f: &mut FunctionBuilder<'t>, x: Value) -> Value {
	let x = f.ins().fdemote(F32, x);
	let exponent = f.ins().iconst(I32, 0x7F800000);
	let mantissa = f.ins().iconst(I32, 0x007FFFFF);
	let i = f.ins()./*raw_?*/bitcast(I32, x);
	let e = f.ins().fcvt_from_sint(I32, f.ins().isub(f.ins().ushr_imm(f.ins().band(i, exponent), 23), f.ins().iconst(I32, 127)));
	let _1 = f.ins().f32const(1.);
	let m = f.ins().bor(f.ins().bitcast(I32, f.ins().band(i, mantissa)), _1);
	let c = [2.28330284476918490682, -1.04913055217340124191, 0.204446009836232697516].map(|c| f.ins().f32const(c));
	let p = f.ins().fma(f.ins().fma(c[2], m, c[1]), m, c[0]);
	let p = f.ins().fmul(p, f.ins().fsub(m, _1)); //?
	f.ins().fpromote(F64, f.ins().fadd(p, e))
}

use std::f64::consts::LN_2;

fn log_arrhenius(f: &mut FunctionBuilder<'t>, RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: RateConstant,
														rcpT: Value, logT: Value) -> Value {
	f.ins().fma(f.ins().f64const(-activation_temperature/LN_2), rcpT, f.ins().fma(f.ins().f64const(temperature_exponent), logT, f.ins().f64const(log_preexponential_factor/LN_2)))
}

fn fdot(f: &mut FunctionBuilder<'t>, iter: impl IntoIterator<Item=(Value, Value)>) -> Value {
	let iter = iter.into_iter();
	let mut sum = {let (a,b) = iter.next().unwrap(); f.ins().fmul(a, b)};
	for (a,b) in iter { sum = f.ins().fma(a, b, sum); }
	sum
}

fn dot<T>(f: &mut FunctionBuilder<'t>, iter: impl IntoIterator<Item=(T, Value)>, mut sum: Option<Value>) -> Option<Value>
where T: num::IsZero + num::IsOne + num::IsMinusOne + Into<f64> {
	for (c,v) in iter.into_iter() {
		//use num::{IsZero, IsOne, IsMinusOne};
		if c.is_zero() {}
		else if c.is_one() { sum = Some(match sum { Some(sum) => f.ins().fadd(sum, v), None => v}); }
		else if c.is_minus_one() { sum = Some(match sum { Some(sum) => f.ins().fsub(sum, v), None => v}); }
		else { let c = f.ins().f64const(c.into()); sum = Some(match sum { Some(sum) => f.ins().fma(c, v, sum), None => f.ins().fmul(c, v) }); }
	}
	sum
}


impl ReactionModel {
fn efficiency(&self, f: &mut FunctionBuilder<'t>, _1: Value, logT: Value, rcpT: Value, mT: Value, mrcpT: Value, concentrations: &[Value], log_k_inf: Value) -> Value {
	use ReactionModel::*; match self {
		Elementary => _1,
		ThreeBody{efficiencies} => { dot(f, efficiencies.iter().copied().zip(concentrations.iter().copied()), None).unwrap() },
		PressureModification{efficiencies, k0} => {
			let Pr = f.ins().fmul(dot(f, efficiencies.iter().copied().zip(concentrations.iter().copied()), None).unwrap(), exp2(f, f.ins().fsub(log_arrhenius(f, *k0, rcpT, logT), log_k_inf)));
			f.ins().fdiv(Pr, f.ins().fadd(_1, Pr))
		}
		Falloff{efficiencies, k0, troe} => {
			let dot = dot(f, efficiencies.iter().copied().zip(concentrations.iter().copied()), None).unwrap();
			let Pr = f.ins().fmul(dot, exp2(f, f.ins().fsub(log_arrhenius(f, *k0, rcpT, logT), log_k_inf)));
			let model::Troe{A, T3, T1, T2} = *troe;
			let logFcent = log2(f, f.ins().fma(f.ins().f64const(1.-A), exp2(f, f.ins().fdiv(mT, f.ins().f64const(T3))), f.ins().fma(f.ins().f64const(A), exp2(f, f.ins().fdiv(mT, f.ins().f64const(T1))), exp2(f, f.ins().fmul(mrcpT, f.ins().f64const(T2))))));
			let C =f.ins().fma(f.ins().f64const(-0.67), logFcent, f.ins().f64const(-0.4*f64::log2(10.)));
			let N = f.ins().fma(f.ins().f64const(-1.27), logFcent, f.ins().f64const(0.75*f64::log2(10.)));
			let logPr𐊛C = f.ins().fadd(log2(f, Pr), C);
			let f1 = f.ins().fdiv(logPr𐊛C, f.ins().fma(f.ins().f64const(-0.14), logPr𐊛C, N));
			let F = exp2(f, f.ins().fdiv(logFcent, f.ins().fma(f1, f1, _1)));
			f.ins().fdiv(Pr, f.ins().fma(Pr, F, F))
		}
	}
}
}

impl Model {
pub fn rate<const CONSTANT: Property>(&self) -> fn(Constant<CONSTANT>, StateVector<CONSTANT>, &mut Derivative<CONSTANT>) {
	let builder = JITBuilder::new(cranelift_module::default_libcall_names());
	let builder_context = FunctionBuilderContext::new();
	let module = JITModule::new(builder);
  let context = module.make_context();
	let PTR = module.target_config().pointer_type();
  let params = [("constant", F64), ("state", PTR), ("rate", PTR)];
	context.func.signature.params = params.iter().map(|(_,r#type)| AbiParam::new(*r#type)).collect();
	let mut builder = FunctionBuilder::new(&mut context.func, &mut builder_context);
	let entry_block = builder.create_block();
	builder.append_block_params_for_function_params(entry_block);
	builder.switch_to_block(entry_block);
	builder.seal_block(entry_block);
	/*let mut variables = std::collections::HashMap::new();
	let declare_variable = {
		let mut index = 0;
		|builder: FunctionBuilder, /*variables: &mut HashMap<String, Variable>,*/ name: &str, r#type: Type| -> Variable {
			let variable = Variable::new(index);
			if !variables.contains_key(name) {
					variables.insert(name.into(), variable);
					builder.declare_var(variable, r#type);
					index += 1;
			}
			variable
		}
	};*/
	//use iter::array_from_iter as from_iter;
	let [constant, state, rate]: [Value; 3] = builder.block_params(entry_block).try_into().unwrap(); /*from_iter(params.into_iter().enumerate().map(|(i, (name, r#type))|{
		let var = declare_variable(builder, name, *r#type); builder.def_var(var, builder.block_params(entry_block)[i]); var
	}));*/
	let flags = MemFlags::new();
	let ref mut f = builder;
	let _1 = f.ins().f64const(1.);
	let _m1 = f.ins().f64const(-1.);
	use std::mem::size_of;
	let T = f.ins().load(F64, flags, state, 0*size_of::<f64>() as i32);
	let logT = log2(f, T);
	let rcpT = f.ins().fdiv(_1, T);
	let T2 = f.ins().fmul(T, T);
	let T3 = f.ins().fmul(T2, T);
	let T4 = f.ins().fmul(T3, T);
	let mT = f.ins().fneg(T);
	let mrcpT = f.ins().fneg(rcpT);
	let Self{species: Species{molar_mass: W, thermodynamics, heat_capacity_ratio, ..}, reactions} = self;
	let a = thermodynamics.iter().map(|s| s.0[1]).collect(): Box<_>;
	use std::array::IntoIter;
	let G_RT =
		a.iter().map(|a| dot(f,
			IntoIter::new([(a[5], rcpT), (a[0]*LN_2, logT), (-a[1]/2., T), ((1./3.-1./2.)*a[2], T2), ((1./4.-1./3.)*a[3], T3), ((1./5.-1./4.)*a[4]/5., T4)]),
			Some(f.ins().f64const(a[0]-a[6]))).unwrap()).collect(): Box<_>;
	let logP0_RT = f.ins().fsub(f.ins().f64const(f64::ln(NASA7::reference_pressure)), logT);
	let variable = f.ins().load(F64, flags, state, 1*size_of::<f64>() as i32);
	let (pressure, volume) = {use Property::*; match CONSTANT {Pressure => (constant, variable), Volume => (variable, constant)}};
	let C = f.ins().fdiv(pressure, T); // n/V = P/kT
	let len = self.len();
	let amounts = (0..len-1).map(|i| f.ins().load(F64, flags, state, ((2+i)*size_of::<f64>()) as i32)).collect(): Box<_>;
	let rcpV = f.ins().fdiv(_1, volume);
	let concentrations = amounts.iter().map(|&n| f.ins().fmul(n/*.max(0.)*/, rcpV)).collect(): Box<_>;
	let Ca = f.ins().fsub(C, dot(f, std::iter::repeat(1.).zip(concentrations.iter().copied()), None).unwrap());
	//if Ca < 0. { dbg!(T, C, concentrations, Ca); throw!(); }
	let ref concentrations = [&concentrations as &[_],&[Ca]].concat();
	let ref log_concentrations = concentrations.iter().map(|&x| log2(f, x)).collect(): Box<_>;
	let mut dtω = (0..len-1).map(|_| None).collect(): Box<_>;
	for Reaction{reactants, products, net, Σnet, rate_constant, model, ..} in reactions.iter() {
		let log_k_inf = log_arrhenius(f, *rate_constant, rcpT, logT);
		let c = model.efficiency(f, _1, logT, rcpT, mT, mrcpT, concentrations, log_k_inf); // todo: CSE
		let Rf = exp2(f, dot(f, reactants.iter().copied().zip(log_concentrations.iter().copied()), Some(log_k_inf)).unwrap());
		let m_log_equilibrium_constant = dot(f, net.iter().chain(IntoIter::new([Σnet])).copied().zip(G_RT.iter().chain(&[logP0_RT]).copied()), None).unwrap();
		let Rr = exp2(f, dot(f, products.iter().copied().zip(log_concentrations.iter().copied()), Some(f.ins().fadd(log_k_inf, m_log_equilibrium_constant))).unwrap());
		//assert!(Rr.is_finite(), "{} {} {}", $dot_products(log_concentrations), log_k_inf, log_equilibrium_constant);
		let R = f.ins().fsub(Rf, Rr);
		let cR = f.ins().fmul(c, R);
		for (index, &ν) in net.iter().enumerate() {
			let ref mut dtω = dtω[index];
			match ν {
					0 => {},
					1 => match dtω { None => *dtω = Some(cR), Some(dtω) => *dtω = f.ins().fadd(*dtω, cR) }
					-1 => match dtω { None => *dtω = Some(cR), Some(dtω) => *dtω = f.ins().fsub(*dtω, cR) }
					ν => match dtω { None => *dtω = Some(f.ins().fmul(f.ins().f64const(ν as f64), cR)), Some(dtω) => *dtω = f.ins().fma(f.ins().f64const(ν as f64), cR, *dtω) }
			}
		}
	}
	let a = {use Property::*; match CONSTANT {Pressure => a, Volume => a.iter().zip(heat_capacity_ratio.iter()).map(|(a,γ)| a.map(|a| a / γ)).collect()}};
	let Cc/*Cp|Cv*/ = a.iter().map(|a| dot(f, IntoIter::new([(a[1], T), (a[2], T2), (a[3], T3), (a[4], T4)]), Some(f.ins().f64const(a[0]))).unwrap());//.collect();
	let m_rcp_ΣCCc = f.ins().fdiv(_m1, fdot(f, concentrations.iter().copied().zip(Cc)));
	let E_T/*H/T|U/T*/ = a.iter().map(|a| dot(f, IntoIter::new([(a[5], rcpT), (a[1]/2., T), (a[2]/3., T2), (a[3]/4., T3), (a[4]/5., T4)]), Some(f.ins().f64const(a[0]))).unwrap());//.collect();
	let dtω = dtω.into_iter().map(|dtω| dtω.unwrap()).collect(): Box<_>;
	let dtT_T = f.ins().fmul(m_rcp_ΣCCc, fdot(f, dtω.iter().copied().zip(E_T)));
	f.ins().store(flags, rate, f.ins().fmul(dtT_T, T), 0*size_of::<f64>() as i32);
	let R_S_Tdtn = f.ins().fmul(f.ins().fdiv(T, pressure), dot(f, W[0..len-1].iter().map(|w| 1. - w/W[len-1]).zip(dtω.iter().copied()), None).unwrap()); // R/A Tdtn (constant pressure: A=V, constant volume: A=P)
	let dtS_S = f.ins().fadd(R_S_Tdtn, dtT_T);
	f.ins().store(flags, rate, f.ins().fmul(dtS_S, variable), 1*size_of::<f64>() as i32);
	let dtn = dtω.into_iter().map(|&dtω| f.ins().fmul(volume, dtω)).collect(): Box<_>;
	for (i, &dtn) in dtn.into_iter().enumerate() { f.ins().store(flags, rate, dtn, ((2+i)*size_of::<f64>()) as i32); }
	builder.finalize();
  let id = module.declare_function(&"", Linkage::Export, &context.func.signature).unwrap();
	module.define_function(id, &mut context, &mut codegen::binemit::NullTrapSink{}).unwrap();
	module.finalize_definitions();
	std::mem::transmute::<_,fn(Constant<CONSTANT>, StateVector<CONSTANT>, &mut Derivative<CONSTANT>)>(module.get_finalized_function(id))
}
}
