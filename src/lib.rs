#![feature(const_generics, const_generics_defaults, const_evaluatable_checked, type_ascription, array_methods, in_band_lifetimes, once_cell, array_map, map_into_keys_values, bindings_after_at, destructuring_assignment, trait_alias, non_ascii_idents)]
#![allow(incomplete_features, non_snake_case,confusable_idents, mixed_script_confusables, non_upper_case_globals, unused_imports, uncommon_codepoints)]
pub use system::default;
mod transport; pub use transport::{TransportPolynomials, Transport};
pub mod reaction; pub use reaction::{Reaction, Model, RateConstant};
use {std::f64::consts::PI as π, num::{sq, cb, sqrt, log, pow, powi}};
use iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoCopied, Enumerate, IntoChain, map}, zip, map, eval, vec::{self, eval, Dot, generate, Scale, Sub}};
use system::{Map, Element, Troe};

pub const K : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)

#[derive(Debug)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure : f64 = 101325. / NA; // 1 atm
	const T_split : f64 = 1000.*K;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn specific_enthalpy(&self, T: f64) -> f64 { let a = self.a(T); a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T } // /R
	pub fn specific_enthalpy_T(&self, T: f64) -> f64 { let a = self.a(T); a[5]/T+a[0]+a[1]/2.*T+a[2]/3.*T*T+a[3]/4.*T*T*T+a[4]/5.*T*T*T*T } // /RT
	pub fn specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); a[6]+a[0]*log(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T } // /R
	//fn dT_specific_heat_capacity(&self, T: f64) -> f64 { 	let a = self.a(T); kB*NA * (a[1]+2.*a[2]*T+3.*a[3]*T*T+4.*a[4]*T*T*T) } // /R
	//fn dT_Gibbs_free_energy(&self, T: f64) -> f64 { let a = self.a(T); (1.-a[0])/T - a[1]/2. - a[2]/12.*T - a[3]/36.*T*T - a[4]/80.*T*T*T - a[5]/(T*T) } // dT((H-TS)/RT)
}

pub struct Species<const S: usize> {
	pub molar_mass: [f64; S],
	pub thermodynamics: [NASA7; S],
	diameter: [f64; S],
	well_depth_J: [f64; S],
	polarizability: [f64; S],
	permanent_dipole_moment: [f64; S],
	rotational_relaxation: [f64; S],
	internal_degrees_of_freedom: [f64; S],
	heat_capacity_ratio: [f64; S],
}

pub struct System<const S: usize> where [(); S-1]: {
	pub species: Species<S>,
	pub reactions: Box<[Reaction<S>]>, // .net[S-1]
	pub transport_polynomials: TransportPolynomials<S>,
}

#[derive(Clone, Copy, Debug, PartialEq)] pub struct State<const S: usize> {
	pub temperature: f64,
	pub pressure: f64,
	pub volume: f64,
	pub amounts: [f64; S]
}

pub struct Simulation<'t, const S: usize> where [(); S-1]: {
	pub species_names: [&'t str; S],
	pub system: System<S>,
	pub time_step: f64,
	pub state: State<S>
}

use std::{convert::TryInto, lazy::SyncLazy};
pub static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, C: 12.011, N: 14.0067, O: 15.999, Ar: 39.95}").unwrap().into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

impl From<system::RateConstant> for reaction::RateConstant {
	fn from(system::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: system::RateConstant) -> Self {
		const J_per_cal: f64 = 4.184;
		Self{log_preexponential_factor: log(preexponential_factor)-temperature_exponent*log(K), temperature_exponent, activation_temperature: activation_energy*J_per_cal/NA}
	}
}

mod parse; use parse::*;
// Compile-time selected system
const SPECIES_LEN: usize = parse(env!("SPECIES_LEN"));

impl<const S: usize /*= SPECIES_LEN*/> Simulation<'t, S> where [(); S-1]: {
	pub fn new(system: &'t [u8]) -> ron::Result<Self> where /*'b: 't,*/ [(); S]: {
		let system::System{species, reactions, time_step, state} = ::ron::de::from_bytes(&system)?;
		let species: Vec<_> = species.into();
		use std::convert::TryInto;
		let species : [_; S] = {let len = species.len(); unwrap::unwrap!(species.try_into(), "Compiled for {} species, got {}", S, len)};
		let ref species = species.map(|(k,mut specie):(_,system::Specie)| {
			for T in specie.thermodynamic.temperature_ranges.iter_mut() { *T *= K; } // K->J
			for piece in specie.thermodynamic.pieces.iter_mut() { for (order, a) in piece[0..6].iter_mut().enumerate() { *a /= f64::powi(K, (1+order) as i32); } } // /K^n->/J^n
			(k, specie)
		});
		let molar_mass = eval(species, |(_,s)| s.composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
		let thermodynamics = eval(species, |(_, system::Specie{thermodynamic: system::NASA7{temperature_ranges, pieces},..})| match temperature_ranges[..] {
			[_,Tsplit,_] if Tsplit == NASA7::T_split => NASA7(pieces[..].try_into().unwrap()),
			[min, max] if min < NASA7::T_split && NASA7::T_split < max => NASA7([pieces[0]; 2]),
			ref ranges => panic!("{:?}", ranges),
		});
		let diameter = eval(species, |(_,s)| s.transport.diameter_Å*1e-10);
		let well_depth_J = eval(species, |(_,s)| s.transport.well_depth_K * K);
		use system::Geometry::*;
		let polarizability = eval(species, |(_,s)| if let Linear{polarizability_Å3,..}|Nonlinear{polarizability_Å3,..} = s.transport.geometry { polarizability_Å3*1e-30 } else { 0. });
		let permanent_dipole_moment = eval(species, |(_,s)| if let Nonlinear{permanent_dipole_moment_Debye,..} = s.transport.geometry { permanent_dipole_moment_Debye*Cm_per_Debye } else { 0. });
		let rotational_relaxation = eval(species, |(_,s)| if let Nonlinear{rotational_relaxation,..} = s.transport.geometry { rotational_relaxation } else { 0. });
		let internal_degrees_of_freedom = eval(species, |(_,s)| match s.transport.geometry { Atom => 0., Linear{..} => 1., Nonlinear{..} => 3./2. });
		let heat_capacity_ratio = eval(species, |(_,s)| {
			let f = match s.transport.geometry { Atom => 3., Linear{..} => 5., Nonlinear{..} => 6. };
			1. + 2./f
		});
		let species_names = eval(species, |(name,_)| *name);
		let species_composition = eval(species, |(_,s)| &s.composition);
		let species = Species{molar_mass, thermodynamics, diameter, well_depth_J, polarizability, permanent_dipole_moment, rotational_relaxation, internal_degrees_of_freedom, heat_capacity_ratio};
		let transport_polynomials = species.transport_polynomials();
		let reactions = (reactions.into():Vec<_>).into_iter().map(|system::Reaction{ref equation, rate_constant, model}| {
			for side in equation { for (specie, _) in side { assert!(species_names.contains(&specie), "{}", specie) } }
			let [reactants, products] = eval(equation, |e| eval(species_names, |s| *e.get(s).unwrap_or(&0) as f64));
			let net = iter::vec::Sub::sub(products.prefix(), reactants.prefix());
			{let mut net_composition = Map::new();
				for (s, ν) in net.enumerate() {
					for (element, &count) in species_composition[s] {
						if !net_composition.contains_key(&element) { net_composition.insert(element, 0); }
						*net_composition.get_mut(&element).unwrap() += ν as i8 * count as i8;
					}
				}
				for (_, &ν) in &net_composition { assert!(ν == 0, "{:?} {:?}", net_composition, equation); }}
			let [Σreactants, Σproducts] = [reactants.iter().sum(), products.iter().sum()];
			let Σnet = Σproducts-Σreactants;
			let from = |efficiencies:Map<_,_>| eval(species_names, |specie| *efficiencies.get(specie).unwrap_or(&1.));
			Reaction{
				reactants, products, net, Σreactants, Σproducts, Σnet,
				rate_constant: rate_constant.into(),
				model: {use {system::Model::*, reaction::Model}; match model {
					Elementary => Model::Elementary,
					ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: from(efficiencies)},
					PressureModification{efficiencies, k0} => Model::PressureModification{efficiencies: from(efficiencies), k0: k0.into()},
					Falloff{efficiencies, k0, troe: Troe{A,T3,T1,T2}} => Model::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe: Troe{A,T3:T3*K,T1:T1*K,T2:T2*K}},
				}},
			}
		}).collect();
		let system = System{species, transport_polynomials, reactions};

		let system::State{temperature, pressure, volume, amount_proportions} = state;
		let pressure = pressure/NA;
		let temperature = temperature*K; // K->J
		let amount = pressure * volume / temperature;
		for (specie,_) in &amount_proportions { assert!(species_names.contains(specie)); }
		let amount_proportions = eval(species_names, |specie| *amount_proportions.get(specie).unwrap_or(&0.));
		let amounts = eval(amount_proportions/*.prefix()*/, |amount_proportion| amount/amount_proportions.iter().sum::<f64>() * amount_proportion);

		Ok(Self{
			species_names,
			system,
			time_step,
			state: State{temperature, pressure, volume, amounts}
		})
	}
}

pub fn new(system: &[u8]) -> ron::Result<Simulation<SPECIES_LEN>> { Simulation::new(&system) }

#[cfg(test)] mod test;

#[derive(PartialEq, Eq)] pub enum Property { Pressure, Volume }

impl<const S: usize> State<S> where [(); S-1]:, [(); 2+S-1]: {
	pub fn new<const CONSTANT: Property>(total_amount: f64, thermodynamic_state_constant: f64, u: &reaction::State<CONSTANT, S>) -> Self {
		let u = u.0;
		let amounts: &[_; S-1] = u.suffix();
		let (pressure, volume) = {use Property::*; match CONSTANT {
			Pressure => (thermodynamic_state_constant, u[1]),
			Volume => (u[1], thermodynamic_state_constant)
		}};
		State{temperature: u[0], pressure, volume, amounts: from_iter(amounts.copied().chain([total_amount - iter::into::Sum::<f64>::sum(amounts)]))}
	}
	//pub fn constant<const CONSTANT: Property>(&Self{pressure, volume, ..}: &Self) -> f64 { // arbitrary_self_types
	pub fn constant<const CONSTANT: Property>(&self) -> f64 { let Self{pressure, volume, ..} = self;
		*{use Property::*; match CONSTANT {Pressure => pressure, Volume => volume}}
	}
}

impl<const S: usize> System<S> where [(); S-1]:, [(); 2+S-1]: {
	pub fn rate<const CONSTANT: Property>(&self, state: &State<S>) -> reaction::Derivative<CONSTANT, S> {
		self.rate_and_jacobian(state.constant::<CONSTANT>(), &state.into()).unwrap().0
	}
}

impl<const S: usize> std::fmt::Display for State<S> {
	fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
		let Self{temperature, pressure, volume, amounts} = self;
		write!(fmt, "T: {}, P: {}, V: {}, n: {:?}", temperature/K, pressure*NA, volume, amounts)
	}
}
