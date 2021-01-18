#![allow(incomplete_features)]
#![feature(const_generics, const_evaluatable_checked, type_ascription, array_methods, non_ascii_idents,in_band_lifetimes,once_cell,array_map,map_into_keys_values,bindings_after_at,destructuring_assignment)]
#![feature(trait_alias)]
#![allow(non_snake_case,confusable_idents,mixed_script_confusables,non_upper_case_globals,unused_imports,uncommon_codepoints)]
pub mod ron;
mod transport; pub use transport::{TransportPolynomials, Transport};
mod reaction; pub use reaction::{Reaction, Model, RateConstant};
use {std::f64::consts::PI as π, num::{sq, cb, sqrt, log, pow, powi}};
use {iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoCopied, Enumerate, IntoChain, map}, zip, map, eval, vec::{self, eval, Dot, generate, Scale, Sub}}, self::ron::{Map, Element, Troe}};

pub const kB : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)

#[derive(Debug)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure_R : f64 = 101325. / (kB*NA); // 1 atm
	const T_split : f64 = 1000.;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn specific_enthalpy(&self, T: f64) -> f64 { let a = self.a(T); a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T } // /R
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
}

pub struct System<const S: usize> where [(); S-1]: {
	pub species: Species<S>,
	pub reactions: Box<[Reaction<S>]>, // .net[S-1]
	pub transport_polynomials: TransportPolynomials<S>,
}
impl<const S: usize> System<S> where [(); S-1]: {
	const volume : f64 = 1.;
}

#[derive(Clone, Copy)] pub struct State<const S: usize> {
	pub temperature: f64,
	pub amounts: [f64; S]
}

pub struct Simulation<'t, const S: usize> where [(); S-1]: {
	pub species_names: [&'t str; S],
	pub system: System<S>,
	pub time_step: f64,
	pub pressure_R: f64,
	pub state: State<S>
}

use std::{convert::TryInto, lazy::SyncLazy};
pub static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	::ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, C: 12.011, O: 15.999, Ar: 39.95}").unwrap().into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

impl From<ron::RateConstant> for reaction::RateConstant { fn from(ron::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: ron::RateConstant) -> Self {
	const J_per_cal: f64 = 4.184;
	Self{log_preexponential_factor: log(preexponential_factor), temperature_exponent, activation_temperature: activation_energy*J_per_cal/(kB*NA)}
}}

impl<const S: usize> Simulation<'t, S> where [(); S-1]: {
	pub fn new(system: &'b [u8]) -> ::ron::Result<Self> where 'b: 't, [(); S]: {
		let ron::System{species, reactions, time_step, state} = ::ron::de::from_bytes(&system)?;
		let species: Vec<_> = species.into();
		use {std::convert::TryInto, iter::into::IntoMap};
		let ref species : [_; S] = {let len = species.len(); unwrap::unwrap!(species.try_into(), "Compiled for {} species, got {}", S, len)};
		let molar_mass = eval(species, |(_,s)| s.composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
		let thermodynamics = eval(species, |(_,ron::Specie{thermodynamic: ron::NASA7{temperature_ranges, pieces},..})| match temperature_ranges[..] {
			[_,Tsplit,_] if Tsplit == NASA7::T_split => NASA7(pieces[..].try_into().unwrap()),
			[min, max] if min < NASA7::T_split && NASA7::T_split < max => NASA7([pieces[0]; 2]),
			ref ranges => panic!("{:?}", ranges),
		});
		let diameter = eval(species, |(_,s)| s.transport.diameter_Å*1e-10);
		let well_depth_J = eval(species, |(_,s)| s.transport.well_depth_K * kB);
		use self::ron::Geometry::*;
		let polarizability = eval(species, |(_,s)| if let Linear{polarizability_Å3,..}|Nonlinear{polarizability_Å3,..} = s.transport.geometry { polarizability_Å3*1e-30 } else { 0. });
		let permanent_dipole_moment = eval(species, |(_,s)| if let Nonlinear{permanent_dipole_moment_Debye,..} = s.transport.geometry { permanent_dipole_moment_Debye*Cm_per_Debye } else { 0. });
		let rotational_relaxation = eval(species, |(_,s)| if let Nonlinear{rotational_relaxation,..} = s.transport.geometry { rotational_relaxation } else { 0. });
		let internal_degrees_of_freedom = eval(species, |(_,s)| match s.transport.geometry { Atom => 0., Linear{..} => 1., Nonlinear{..} => 3./2. });
		let species_names = eval(species, |(name,_)| *name);
		let species = Species{molar_mass, thermodynamics, diameter, well_depth_J, polarizability, permanent_dipole_moment, rotational_relaxation, internal_degrees_of_freedom};
		let transport_polynomials = species.transport_polynomials();
		let reactions = iter::into::Collect::collect(reactions.map(|self::ron::Reaction{ref equation, rate_constant, model}| {
			for side in equation { for (specie, _) in side { assert!(species_names.contains(&specie), "{}", specie) } }
			let [reactants, products] = eval(equation, |e| eval(species_names, |s| *e.get(s).unwrap_or(&0) as f64));
			let net = iter::vec::Sub::sub(products.prefix(), reactants.prefix());
			let [Σreactants, Σproducts] = [reactants.iter().sum(), products.iter().sum()];
			let Σnet = Σproducts-Σreactants;
			let from = |efficiencies:Map<_,_>| eval(species_names, |specie| *efficiencies.get(specie).unwrap_or(&1.));
			Reaction{
				reactants, products, net, Σreactants, Σproducts, Σnet,
				rate_constant: rate_constant.into(),
				model: {use self::{ron::Model::*, reaction::Model}; match model {
					Elementary => Model::Elementary,
					ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: from(efficiencies)},
					PressureModification{efficiencies, k0} => Model::PressureModification{efficiencies: from(efficiencies), k0: k0.into()},
					Falloff{efficiencies, k0, troe} => Model::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe},
				}},
			}
		}));
		let system = System{species, transport_polynomials, reactions};

		let ron::State{temperature, pressure, mole_proportions} = state;
		let pressure_R = pressure / (kB*NA);
		let amount = pressure_R / temperature * System::<S>::volume;
		let mole_proportions = eval(species_names, |specie| *mole_proportions.get(specie).unwrap_or(&0.));
		let amounts = eval(mole_proportions/*.prefix()*/, |mole_proportion| amount/mole_proportions.iter().sum::<f64>() * mole_proportion);

		Ok(Self{
			species_names,
			system,
			time_step, pressure_R,
			state: State{temperature, amounts}
		})
	}
}

//internal compiler error: compiler/rustc_middle/src/ich/impls_ty.rs:94:17: StableHasher: unexpected region '_#0r
/*impl<const S: usize> From<&State<S>> for [f64; 1+S-1] where [(); S-1]: { fn from(State{temperature, amounts}: &State<S>) -> Self {
	use iter::{array_from_iter as from_iter, into::IntoChain}; from_iter([*temperature].chain(amounts[..S-1].try_into().unwrap():[_;S-1]))
} }*/
impl<const S: usize> From<State<S>> for [f64; 1+S-1] where [(); S-1]: { fn from(State{temperature, amounts}: State<S>) -> Self {
	use iter::{array_from_iter as from_iter, into	::IntoChain}; from_iter([temperature].chain(amounts[..S-1].try_into().unwrap():[_;S-1]))
} }
impl<const S: usize> State<S> {
	pub fn new(total_amount: f64, u: [f64; 1+S-1]) -> Self where [(); S-1]: {
		let amounts: &[_; S-1] = u.suffix();
		State{temperature: u[0], amounts: from_iter(amounts.copied().chain([total_amount - iter::into::Sum::<f64>::sum(amounts)]))}
	}
}

pub trait Error { fn error(&self, o: &Self) -> f64; }
impl Error for f64 { fn error(&self, o: &Self) -> f64 { f64::abs(self-o) } }
impl<const N: usize> Error for [f64; N] { fn error(&self, o: &Self) -> f64 { self.iter().zip(o).map(|(s,o)| s.error(o)).max_by(|s,o| s.partial_cmp(o).unwrap_or_else(||panic!("{} {}", s, o))).unwrap() } }
impl<const S: usize> Error for Transport<S> {
		fn error(&self, o: &Self) -> f64 {
		self.viscosity.error(&o.viscosity).max(
		self.thermal_conductivity.error(&o.thermal_conductivity).max(
		self.mixture_averaged_thermal_diffusion_coefficients.error(&o.mixture_averaged_thermal_diffusion_coefficients) ))
	}
}
