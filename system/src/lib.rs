#![feature(const_generics, const_evaluatable_checked, non_ascii_idents, type_ascription, once_cell)]
//, type_ascription, array_methods, in_band_lifetimes, , array_map, map_into_keys_values, bindings_after_at, destructuring_assignment, trait_alias)]
//#![allow(non_snake_case, , mixed_script_confusables, non_upper_case_globals, unused_imports, uncommon_codepoints)]
#![allow(incomplete_features, non_upper_case_globals, non_snake_case, confusable_idents)]
use std::ops::Deref;
//use {std::f64::consts::PI as π, num::{sq, cb, sqrt, log, pow, powi}};
use num::log;
//use iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoCopied, Enumerate, IntoChain, map}, zip, map, eval, vec::{self, eval, Dot, generate, Scale, Sub}};
//use model::{/*Map, Element,*/ Troe};
//mod transport; pub use transport::{TransportPolynomials, Transport};
//pub mod reaction; pub use reaction::{Reaction, Model, RateConstant};*/

pub const K : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;

#[derive(PartialEq, Debug, /*Eq*/)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure : f64 = 101325. / NA; // 1 atm
	pub const T_split : f64 = 1000.*K;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn specific_enthalpy(&self, T: f64) -> f64 { let a = self.a(T); a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T } // /R
	pub fn specific_enthalpy_T(&self, T: f64) -> f64 { let a = self.a(T); a[5]/T+a[0]+a[1]/2.*T+a[2]/3.*T*T+a[3]/4.*T*T*T+a[4]/5.*T*T*T*T } // /RT
	pub fn specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); a[6]+a[0]*log(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T } // /R
	//fn dT_specific_heat_capacity(&self, T: f64) -> f64 { 	let a = self.a(T); kB*NA * (a[1]+2.*a[2]*T+3.*a[3]*T*T+4.*a[4]*T*T*T) } // /R
	//fn dT_Gibbs_free_energy(&self, T: f64) -> f64 { let a = self.a(T); (1.-a[0])/T - a[1]/2. - a[2]/12.*T - a[3]/36.*T*T - a[4]/80.*T*T*T - a[5]/(T*T) } // dT((H-TS)/RT)
}

#[derive(PartialEq, /*Eq,*/ Debug, Clone, Copy)] pub struct RateConstant {
	pub log_preexponential_factor: f64,
	pub temperature_exponent: f64,
	pub activation_temperature: f64
}

pub fn log_arrhenius(RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: RateConstant, T: f64) -> f64 {
	log_preexponential_factor + temperature_exponent*num::log(T) - activation_temperature*(1./T)
}

pub use model;

impl From<model::RateConstant> for RateConstant {
	fn from(model::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: model::RateConstant) -> Self {
		const J_per_cal: f64 = 4.184;
		Self{log_preexponential_factor: log(preexponential_factor)-temperature_exponent*log(K), temperature_exponent, activation_temperature: activation_energy*J_per_cal/NA}
	}
}

#[derive(PartialEq, Eq)] pub enum Property { Pressure, Volume }

#[derive(Debug)] pub struct State<const CONSTANT: Property, const S: usize>(pub [f64; 2+S-1]) where [(); 2+S-1]:;
pub type Derivative<const CONSTANT: Property, const S: usize> = State<CONSTANT, S>;

impl<const CONSTANT: Property, const S: usize> std::ops::Deref for State<CONSTANT, S> where [(); 2+S-1]: {
	type Target = [f64; 2+S-1]; fn deref(&self) -> &Self::Target { &self.0 }
}

use linear_map::LinearMap as Map;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)
use model::{/*Map,*/ Element, Troe};

use std::{convert::TryInto, lazy::SyncLazy};
static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, C: 12.011, N: 14.0067, O: 15.999, Ar: 39.95}").unwrap()
	.into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

#[derive(PartialEq/*, Eq*/)] #[allow(dead_code)] pub struct Species<const S: usize> {
	pub molar_mass: [f64; S],
	pub thermodynamics: [NASA7; S],
	pub diameter: [f64; S],
	pub well_depth_J: [f64; S],
	pub polarizability: [f64; S],
	pub permanent_dipole_moment: [f64; S],
	pub rotational_relaxation: [f64; S],
	pub internal_degrees_of_freedom: [f64; S],
	pub heat_capacity_ratio: [f64; S],
}

#[derive(Debug)] pub struct Species_ {
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

#[derive(Debug)] pub enum Model {
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
	pub model: Model,
}

#[derive(Debug)] pub struct System {
	pub species: Species_,
	pub reactions: Box<[Reaction]>,
	//pub transport_polynomials: TransportPolynomials_,
}

impl System {
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
				Elementary => Model::Elementary,
				ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: from(efficiencies)},
				PressureModification{efficiencies, k0} => Model::PressureModification{efficiencies: from(efficiencies), k0: k0.into()},
				Falloff{efficiencies, k0, troe: Troe{A,T3,T1,T2}} => Model::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe: Troe{A,T3:T3*K,T1:T1*K,T2:T2*K}},
			}},
		}
	});
	System{species, /*transport_polynomials,*/ reactions}
}
pub fn len(&self) -> usize { self.species.molar_mass.len() }
}
