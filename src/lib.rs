#![feature(nll, const_generics, const_evaluatable_checked, once_cell, in_band_lifetimes, array_map, trait_alias, unboxed_closures, fn_traits, array_methods, bindings_after_at, associated_type_bounds, default_free_fn)]
#![allow(incomplete_features, non_upper_case_globals, non_snake_case, confusable_idents, uncommon_codepoints)]
#![allow(unused_variables, dead_code)]

pub const K : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)

pub mod model;
use model::Element;

#[derive(PartialEq)] pub struct NASA7 {
	pub temperature_split : f64,
	pub pieces: [[f64; 7]; 2],
}

impl NASA7 {
	pub const reference_pressure : f64 = 101325. / (K*NA);
	pub fn piece(&self, T: f64) -> &[f64; 7] { &self.pieces[if T < self.temperature_split { 0 } else { 1 }] }
	pub fn molar_heat_capacity_at_constant_pressure_R(&self, T: f64) -> f64 { let a = self.piece(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
}

use {std::{convert::TryInto, lazy::SyncLazy}, linear_map::LinearMap as Map};
static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, C: 12.011, N: 14.0067, O: 15.999, Ar: 39.95}").unwrap()
	.into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

pub struct Species {
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
impl Species {
	pub fn new(species: &Map<&'t str, model::Specie>) -> (Box<[&'t str]>, Self) {
		/*let species = eval(species, |(k,specie):&(_,model::Specie)| {
			let mut specie = specie.clone();
			for T in specie.thermodynamic.temperature_ranges.iter_mut() { *T *= K; } // K->J
			for piece in specie.thermodynamic.pieces.iter_mut() { for (order, a) in piece[0..6].iter_mut().enumerate() { *a /= f64::powi(K, (1+order) as i32); } } // /K^n->/J^n
			(k, specie)
		});
		use std::ops::Deref;
		let species = species.deref();*/
		use iter::map;
		let molar_mass = map(species, |(_,s)| s.composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
		let thermodynamics = map(species, |(_, model::Specie{thermodynamic: model::NASA7{temperature_ranges, pieces},..})| match temperature_ranges[..] {
			[_,temperature_split,_] => NASA7{temperature_split, pieces: pieces[..].try_into().unwrap()},
			[min, max] => NASA7{temperature_split: f64::INFINITY, pieces: [pieces[0]; 2]},
			ref ranges => panic!("{:?}", ranges),
		});
		let diameter = map(species, |(_,s)| s.transport.diameter_Å*1e-10);
		let well_depth_J = map(species, |(_,s)| s.transport.well_depth_K * K);
		use model::Geometry::*;
		let polarizability = map(species, |(_,s)| if let Linear{polarizability_Å3,..}|Nonlinear{polarizability_Å3,..} = s.transport.geometry { polarizability_Å3*1e-30 } else { 0. });
		let permanent_dipole_moment = map(species, |(_,s)|
			if let Nonlinear{permanent_dipole_moment_Debye,..} = s.transport.geometry { permanent_dipole_moment_Debye*Cm_per_Debye } else { 0. });
		let rotational_relaxation = map(species, |(_,s)| if let Nonlinear{rotational_relaxation,..} = s.transport.geometry { rotational_relaxation } else { 0. });
		let internal_degrees_of_freedom = map(species, |(_,s)| match s.transport.geometry { Atom => 0., Linear{..} => 1., Nonlinear{..} => 3./2. });
		let heat_capacity_ratio = map(species, |(_,s)| {
			let f = match s.transport.geometry { Atom => 3., Linear{..} => 5., Nonlinear{..} => 6. };
			1. + 2./f
		});
		let species_names = map(species, |(name,_)| *name);
		let species_composition = map(species, |(_,s)| &s.composition);
		(species_names,
			Species{molar_mass, thermodynamics, diameter, well_depth_J, polarizability, permanent_dipole_moment, rotational_relaxation, internal_degrees_of_freedom,
										heat_capacity_ratio})
	}
	pub fn len(&self) -> usize { self.molar_mass.len() }
}

pub struct State {
    pub temperature: f64,
    pub pressure_R: f64,
    pub volume: f64,
    pub amounts: Box<[f64]>
}

pub fn initial_state(model::Model{species, state, time_step, ..}: &model::Model<'t>) -> State {
	let species_names = iter::map(species, |(name,_)| *name);

	let model::State{temperature, pressure, volume, amount_proportions} = state;
	let pressure_R = pressure/(K*NA);
	let temperature = *temperature; //K*: K->J
	let amount = pressure_R * volume / temperature;
	for (specie,_) in amount_proportions { assert!(species_names.contains(specie)); }
	let amount_proportions = iter::map(&*species_names, |specie| *amount_proportions.get(specie).unwrap_or(&0.));
	let amounts = amount_proportions.iter().map(|amount_proportion| amount * amount_proportion/amount_proportions.iter().sum::<f64>()).collect();

	State{temperature, pressure_R, volume: *volume, amounts}
}

#[cfg(feature= "transport")] pub mod transport;
#[cfg(feature= "reaction")] pub mod reaction;
