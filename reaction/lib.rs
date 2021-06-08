#![feature(trait_alias)]#![allow(non_snake_case)]
use {fehler::throws, anyhow::Error};

use combustion::*;
pub struct Simulation<'t> {
	pub species_names: Box<[&'t str]>,
	pub function: Function,
	pub states: Box<[f64]>,
	pub rates: Box<[f64]>,
	pub pressure_Pa_R: f64,
	pub temperature: f64,
	pub mass_rate: f64,
	pub energy_rate_R: f64
}

impl<'t> Simulation<'t> {
#[throws] pub fn new(model: &'t [u8], states_len: usize) -> Self {
	let model = model::Model::new(&model)?;
	let (species_names, species) = combustion::Species::new(&model.species);
	use reaction::*;
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(&species_names, r));
	let function = rate::<_,{Property::Pressure}>(&species, &*reactions, states_len);

	let ref reference_state = initial_state(&model);
	let length = 1f64;
	let velocity = 1f64;
	let time = length / velocity;
	let temperature = reference_state.temperature;
	let total_concentration = reference_state.pressure_R / temperature;
	let total_amount = reference_state.amounts.iter().sum::<f64>();
	let mole_fractions = iter::map(&*reference_state.amounts, |n| n / total_amount);
	let dot = |a:&[_],b:&[_]| iter::dot(a.iter().copied().zip(b.iter().copied()));
	let molar_mass = dot(&*mole_fractions, &*species.molar_mass);
	let density = total_concentration * molar_mass;
	let molar_heat_capacity_R = iter::dot(mole_fractions.iter().copied().zip(
		species.thermodynamics.iter().map(|specie| specie.molar_heat_capacity_at_constant_pressure_R(temperature))
	));
	let mass_rate = density / time ;
	let energy_rate_R = (molar_heat_capacity_R * reference_state.pressure_R) / time;

	let ref state = initial_state(&model);
	let pressure_Pa_R = state.pressure_R;
	let total_amount = state.amounts.iter().sum::<f64>();
	let mole_fractions = iter::map(&*state.amounts, |n| n / total_amount);
	let molar_mass = dot(&*mole_fractions, &*species.molar_mass);
	let mass_fractions = iter::map(mole_fractions.iter().zip(species.molar_mass.iter()), |(x,m)| x * m / molar_mass);
	for (&mass_fraction, name) in mass_fractions.iter().zip(species_names.iter()) { if mass_fraction != 0. { println!("{:3} {}",name, mass_fraction) } }
	let state = [&[dbg!(state.temperature / temperature)] as &[_], &*mass_fractions].concat();

	let states = iter::box_collect(state.iter().map(|&s| std::iter::repeat(s as f64).take(states_len)).flatten());
	let rates = vec![f64::NAN; (1/*2*/+species.len()-1)*states_len].into_boxed_slice();
	Self{species_names, function, states, rates, pressure_Pa_R, temperature, mass_rate, energy_rate_R}
}
	pub fn states_len(&self) -> usize { self.rates.len()/(1+self.species_names.len()-1) }
}

pub fn report(species_names: &[&str], rates: &[f64]) {
	#[track_caller] fn all_same(slice: &[f64], stride: usize) -> Box<[f64]> {
		assert_eq!(slice.len()%stride, 0);
		iter::eval(slice.len()/stride, |i| {
			let slice = &slice[i*stride..(i+1)*stride];
			for &v in slice.iter() { assert_eq!(v.to_bits(), slice[0].to_bits()); }
			slice[0]
		})
	}
	let states_len = rates.len()/(1+species_names.len()-1);
	let rates = all_same(&rates, states_len);
	for (&mass_production_rate, name) in rates[1..].iter().zip(species_names.iter()) {
		if mass_production_rate != 0. { println!("{:8} {:e}", name, mass_production_rate); }
	}
	println!("{:8} {:e}", "HRR/HRR0", rates[0]);
}
