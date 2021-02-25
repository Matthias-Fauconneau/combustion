//#![feature(const_generics, const_generics_defaults, const_evaluatable_checked, array_methods, in_band_lifetimes, once_cell, map_into_keys_values, bindings_after_at, destructuring_assignment, trait_alias)]
//#![feature(non_ascii_idents, type_ascription, once_cell, array_map)]
#![feature(const_generics, type_ascription, non_ascii_idents, in_band_lifetimes, const_evaluatable_checked)]
//#![allow(incomplete_features, mixed_script_confusables, unused_imports, uncommon_codepoints)]
#![allow(incomplete_features, confusable_idents, non_snake_case, non_upper_case_globals)]

use std::{convert::TryInto, ops::Deref};
//use {std::f64::consts::PI as π, num::{sq, cb, sqrt, log, pow, powi}};
//use iter::{Prefix, Suffix, array_from_iter as from_iter, into::{IntoCopied, Enumerate, IntoChain, map}, zip, map, eval, vec::{self, eval, Dot, generate, Scale, Sub}};
use iter::{array_from_iter as from_iter, vec::eval, Suffix, into::{IntoCopied, IntoChain}};
/*use std::ops::Deref;
use linear_map::LinearMap as Map;
use system::K;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)
use model::{/*Map,*/ Element, Troe};*/
use system::{NA, K, Property};
//mod transport; pub use transport::{TransportPolynomials, Transport};*/
//pub mod reaction; //pub use reaction::{/*Reaction,*/ Property};
//use system::{RateConstant, NASA7};
//use system::*;

/*#[derive(PartialEq/*, Eq*/)] pub struct System<const S: usize, const R: usize> where [(); S-1]: {
	pub species: Species<S>,
	pub reactions: [Reaction<S>; R],
	//pub transport_polynomials: TransportPolynomials<S>,
}*/

#[derive(Clone, Copy, Debug, PartialEq)] pub struct State<const S: usize> {
	pub temperature: f64,
	pub pressure: f64,
	pub volume: f64,
	pub amounts: [f64; S]
}

impl<const CONSTANT: Property, const S: usize> From<&State<S>> for system::State<CONSTANT, S> where [(); S-1]:, [(); 2+S-1]: {
	fn from(State{temperature, pressure, volume, amounts}: &State<S>) -> Self {
		Self(from_iter([*temperature, *{use Property::*; match CONSTANT {Pressure => volume, Volume => pressure}}].chain(amounts[..S-1].try_into().unwrap():[_;S-1])))
	}
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

impl<const S: usize> State<S> where [(); S-1]:, [(); 2+S-1]: {
	pub fn new<const CONSTANT: Property>(total_amount: f64, thermodynamic_state_constant: f64, u: &system::State<CONSTANT, S>) -> Self {
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

/*impl<const S: usize, const R: usize> System<S, R> where [(); S-1]:, [(); 2+S-1]: {
	pub fn rate<const CONSTANT: Property>(/*const*/ &self, state: &State<S>) -> system::Derivative<CONSTANT, S> {
		self.rate_and_jacobian::<CONSTANT>(state.constant::<CONSTANT>(), &state.into()).unwrap().0
	}
}*/

impl<const S: usize> std::fmt::Display for State<S> {
	fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
		let Self{temperature, pressure, volume, amounts} = self;
		write!(fmt, "T: {}, P: {}, V: {}, n: {:?}", temperature/K, pressure*NA, volume, amounts)
	}
}
