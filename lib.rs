#![feature(once_cell,in_band_lifetimes,array_methods,format_args_capture,associated_type_bounds,trait_alias,default_free_fn,type_ascription,array_zip,unboxed_closures,fn_traits)]
#![allow(non_upper_case_globals,non_snake_case,uncommon_codepoints)]
//#![recursion_limit="9"]
pub mod model;
pub use model::{kB, NA};
const light_speed : f64 = 299_792_458.;
const Cm_per_Debye : f64 = 1e-21 / light_speed; //C·m (Coulomb=A⋅s)

#[derive(PartialEq, Debug)] pub struct NASA7 {
	pub temperature_split : f64,
	pub pieces: [[f64; 7]; 2],
}

impl NASA7 {
	pub const reference_pressure : f64 = 101325. / (kB*NA);
	pub fn piece(&self, T: f64) -> &[f64; 7] { &self.pieces[if T <= self.temperature_split { 0 } else { 1 }] }
	pub fn molar_heat_capacity_at_constant_pressure_R(&self, T: f64) -> f64 { let a = self.piece(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn enthalpy_RT(&self, T: f64) -> f64 { let a = self.piece(T); a[0]+a[1]/2.*T + a[2]/3.*T*T + a[3]/4.*T*T*T + a[4]/5.*T*T*T*T + a[5]/T } // /RT
	pub fn gibbs_RT(&self, T: f64) -> f64 { let a = self.piece(T); a[0]-a[6]-a[0]*f64::ln(T)-a[1]/2.*T+(1./3.-1./2.)*a[2]*T*T+(1./4.-1./3.)*a[3]*T*T*T+(1./5.-1./4.)*a[4]*T*T*T*T+a[5]/T }
}

use {std::lazy::SyncLazy, linear_map::LinearMap as Map, model::Element};
static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(||
	{use Element::*; [(H, 1.008), (He, 4.002602), (C, 12.011), (N, 14.007), (O, 15.999), (F, 18.998403163), (Cl, 35.45), (Ar, 39.95)]}.map(|(e,g)| (e, g/1e3/*kg/g*/)).into_iter().collect()
);

#[derive(Debug)] pub struct Species {
	pub molar_mass: Box<[f64]>,
	pub thermodynamics: Box<[NASA7]>,
	pub diameter: Box<[f64]>,
	pub well_depth_J: Box<[f64]>,
	pub polarizability: Box<[f64]>,
	pub permanent_dipole_moment: Box<[f64]>,
	pub rotational_relaxation: Box<[f64]>,
	pub internal_degrees_of_freedom: Box<[f64]>,
	pub heat_capacity_ratio: Box<[f64]>,
}

use iter::map;

impl Species {
	pub fn new(species: &[(&'t str, model::Specie)]) -> Self {
		let molar_mass = map(species, |(_,s)| s.composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
		let thermodynamics = map(species, |(_, model::Specie{thermodynamic: model::NASA7{temperature_ranges, pieces},..})| match temperature_ranges[..] {
			[_, temperature_split, _] => NASA7{temperature_split, pieces: pieces[..].try_into().unwrap()},
			[_, _] => NASA7{temperature_split: f64::NAN, pieces: [pieces[0]; 2]},
			ref ranges => panic!("{ranges:?}, {species:?}"),
		});
		let diameter = map(species, |(_,s)| s.transport.diameter_Å*1e-10);
		let well_depth_J = map(species, |(_,s)| s.transport.well_depth_K * kB);
		use model::Geometry::*;
		let polarizability = map(species, |(_,s)| if let Linear{polarizability_Å3,..}|Nonlinear{polarizability_Å3,..} = s.transport.geometry { polarizability_Å3*1e-30 } else { 0. });
		let permanent_dipole_moment = map(species, |(_,s)|
			if let Nonlinear{permanent_dipole_moment_Debye,..} = s.transport.geometry { permanent_dipole_moment_Debye*Cm_per_Debye } else { 0. });
		let rotational_relaxation =
			map(species, |(_,s)| if let Linear{rotational_relaxation,..}|Nonlinear{rotational_relaxation,..} = s.transport.geometry { rotational_relaxation } else { 0. });
		let internal_degrees_of_freedom = map(species, |(_,s)| match s.transport.geometry { Atom => 0., Linear{..} => 1., Nonlinear{..} => 3./2. });
		let heat_capacity_ratio = map(species, |(_,s)| 1. + 2. / match s.transport.geometry { Atom => 3., Linear{..} => 5., Nonlinear{..} => 6. });
		Species{molar_mass, thermodynamics, diameter, well_depth_J, polarizability, permanent_dipole_moment, rotational_relaxation, internal_degrees_of_freedom, heat_capacity_ratio}
	}
	pub fn len(&self) -> usize { self.molar_mass.len() }
}

pub struct State {
    pub temperature: f64,
    pub pressure_R: f64,
    pub volume: f64,
    pub amounts: Box<[f64]>
}

pub fn initial_state(model::Model{species, state, ..}: &model::Model<'t>) -> State {
	let model::State{temperature, pressure, volume, amount_proportions} = state;
	let species_names = map(&**species, |(name,_)| *name);
	for (specie,_) in &**amount_proportions { assert!(species_names.contains(specie)); }
	let pressure_R = pressure/(kB*NA);
	let temperature = *temperature;
	let amount = pressure_R * volume / temperature;
	let amount_proportions = map(&*species_names, |specie| *amount_proportions.iter().find(|(s,_)| s==specie).map(|(_,s)| s).unwrap_or(&0.));
	let amounts = map(&*amount_proportions, |amount_proportion| amount * amount_proportion/amount_proportions.iter().sum::<f64>());
	State{temperature, pressure_R, volume: *volume, amounts}
}

pub use model::{RateConstant, Troe};

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

impl Reaction {
	pub fn new(species_names: &[&str], active: usize, model::Reaction{equation, rate_constant, model}: &model::Reaction) -> Self {
		for side in equation { for (specie, _) in side { assert!(species_names.contains(&specie), "{}", specie) } }
		let [reactants, products] = equation.each_ref().map(|e| species_names.iter().map(|&s| *e.get(s).unwrap_or(&0)).collect::<Box<_>>());
		let net = products.into_iter().zip(reactants.into_iter()).take(active).map(|(&a, &b)| a as i8 - b as i8).collect();
		let [Σreactants, Σproducts] = [reactants.iter().sum(), products.iter().sum()];
		let Σnet = Σproducts as i8 - Σreactants as i8;
		let from = |efficiencies:&Map<_,_>| map(species_names, |&specie| *efficiencies.get(specie).unwrap_or(&1.));
		Reaction{
			reactants, products, net, Σreactants, Σproducts, Σnet,
			rate_constant: *rate_constant,
			model: {use model::ReactionModel::*; match model {
				Elementary => ReactionModel::Elementary,
				Irreversible => ReactionModel::Irreversible,
				ThreeBody{efficiencies} => ReactionModel::ThreeBody{efficiencies: from(efficiencies)},
				PressureModification{efficiencies, k0} => ReactionModel::PressureModification{efficiencies: from(efficiencies), k0: *k0},
				Falloff{efficiencies, k0, troe} => ReactionModel::Falloff{efficiencies: from(efficiencies), k0: *k0, troe: *troe},
			}}
		}
	}
}

pub fn new(model: &'t model::Model) -> (Box<[&'t str]>, Species, usize, Box<[Reaction]>, State) {
	let species_names = map(&*model.species, |(name,_)| *name);
	let species = Species::new(&model.species);
	let active = {
		let ref mut iter = species_names.iter().map(
			|specie| model.reactions.iter().any(|model::Reaction{equation,..}| equation[0].get(specie).unwrap_or(&0) != equation[1].get(specie).unwrap_or(&0)) );
		let active = iter.take_while(|is_active| *is_active).count();
		assert!(iter.all(|is_active| !is_active));
		active
	};
	let reactions = map(&*model.reactions, |r| Reaction::new(&species_names, active, r));
	let initial_state = initial_state(&model);
	(species_names, species, active, reactions, initial_state)
}

pub mod reaction;
#[cfg(feature="transport")] pub mod transport;
