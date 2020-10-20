#![feature(in_band_lifetimes, box_syntax)]
#![allow(non_camel_case_types, non_snake_case,non_upper_case_globals)]
mod plot; use plot::real; use plot::Plot;

fn scale(s: real, v: impl Iterator<Item=real>+'t) -> impl Iterator<Item=real>+'t { v.map(move |v| s*v) }
fn recip(x: &'t [real]) -> impl Iterator<Item=real>+'t { x.iter().map(|&x| real::recip(x)) }
fn mul(a: impl Iterator<Item=real>+'t, b: impl Iterator<Item=real>+'t) -> impl Iterator<Item=real>+'t { a.zip(b).map(|(a,b)| a*b) }
fn dot(a: &[real], b: &[real]) -> real { a.iter().zip(b.iter()).map(|(&a,&b)| a*b).sum() }
fn acc(a: &mut [real], b: impl Iterator<Item=real>) { for (a,b) in a.iter_mut().zip(b) { *a += b; } }

use std::collections::BTreeMap as Map;
use serde::Deserialize;
#[derive(Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)] enum Element { H, O, Ar }
#[derive(Deserialize, Debug)] struct InitialState<'t> { temperature: real, pressure: real, volume: real, #[serde(borrow)] mole_proportions: Map<&'t str, real>}
#[derive(Deserialize, Debug)] enum Phase<'t> {
	IdealGas {
		elements: Box<[Element]>,
		species: Box<[&'t str]>,
		#[serde(borrow)] state: InitialState<'t>,
	}
}
#[derive(Deserialize, Debug)] struct NASA7 {
	temperature_ranges: Box<[real]>,
	coefficients: Box<[Box<[real]>]>,
}

#[derive(Deserialize, Debug)] enum Transport {
	Atom { well_depth: real, diameter: real},
	Linear { well_depth: real, diameter: real, polarizability: real, rotational_relaxation: real},
	Nonlinear { well_depth: real, diameter: real, rotational_relaxation: real},
}
#[derive(Deserialize, Debug)] struct Specie {
	composition: Map<Element, u8>,
	specific_heat: NASA7,
	transport: Transport
}

#[derive(Deserialize, Debug)] struct RateConstant {
	#[serde(rename="A")] preexponential_factor: real,
	#[serde(rename="b")] temperature_exponent: real,
	#[serde(rename="Ea")] activation_energy: real
}

#[derive(Deserialize, Debug)] struct Troe { A: real, T3: real, T1: real, T2: real }
#[derive(Deserialize, Debug)] enum Model<'t> {
	Elementary { rate_constants: Box<[RateConstant]> },
	ThreeBody { rate_constants: Box<[RateConstant]>, #[serde(borrow)] efficiencies: Map<&'t str, real> },
	Falloff { k0: RateConstant, k_inf: RateConstant, #[serde(borrow)] efficiencies: Map<&'t str, real>, troe: Troe },
}

#[derive(Deserialize, Debug)] struct Reaction<'t> {
	#[serde(borrow)] equation: (Map<&'t str, u8>, Map<&'t str, u8>),
	model: Model<'t>,
}

#[derive(Deserialize, Debug)] struct Mechanism<'t> {
	time_step: real,
	#[serde(borrow)] phases: Box<[Phase<'t>]>,
	#[serde(borrow)] species: Map<&'t str, Specie>,
	#[serde(borrow)] reactions: Box<[Reaction<'t>]>,
}

const ideal_gas_constant : real = real(8.31446261815324);

impl NASA7 {
	fn at(&self, T: real) -> real {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T)
	}
}

fn arrhenius(&RateConstant{preexponential_factor, temperature_exponent, activation_energy}: &RateConstant, temperature: real) -> real {
	preexponential_factor*temperature.pow(temperature_exponent)*real::exp(-activation_energy/(ideal_gas_constant/real(4.184)*temperature))
}

#[derive(Debug)] enum ModelParameters {
	Elementary { rate_constants: Box<[RateConstant]> },
	ThreeBody { rate_constants: Box<[RateConstant]>, efficiencies: Box<[real]> },
	Falloff { k0: RateConstant, k_inf: RateConstant, efficiencies: Box<[real]>, troe: Troe },
}

fn rate(reactants: (&[usize], &[u8]), model: &ModelParameters, T: real, concentrations: &[real]) -> real {
	let reactants : real = reactants.0.iter().zip(reactants.1.iter()).map(|(&specie, &coefficient)| concentrations[specie].powi(coefficient as i32)).sum();
	use ModelParameters::*;
	reactants * match model {
		Elementary{rate_constants} => rate_constants.iter().map(|rate_constant| arrhenius(rate_constant, T)).sum::<real>(),
		ThreeBody{rate_constants, efficiencies} => rate_constants.iter().map(|rate_constant| arrhenius(rate_constant, T) * reactants).sum::<real>() * dot(efficiencies, concentrations),
		Falloff{k_inf, k0, efficiencies, troe: Troe{A, T3, T1, T2}} => {
			let k_inf = arrhenius(k_inf, T);
			let Pr = arrhenius(k0, T) / k_inf * dot(efficiencies, concentrations);
			let Fcent = (real(1.)-A)*real::exp(-T/T3)+A*real::exp(-T/T1)+real::exp(-T2/T);
			let log10Fcent = real::log10(Fcent);
			let C = real(-0.4)-real(0.67)*log10Fcent;
			let N = real(0.75)-real(1.27)*log10Fcent;
			let f1 = (real::log10(Pr) + C)/(N-real(0.14)*(real::log10(Pr) + C));
			let F = real::exp(real::ln(real(10.))*log10Fcent/(real(1.)+f1*f1));
			k_inf * Pr / (real(1.)+Pr) * F
		}
	}
}

#[fehler::throws(anyhow::Error)] fn main() {
	let standard_atomic_weights : Map<Element, real> = ron::de::from_str("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}")?;
	let mechanism = std::fs::read("H2+O2.ron")?;
	let Mechanism{species, reactions, phases, time_step} = ron::de::from_bytes(&mechanism)?;
	fn from_iter<T>(iter: impl IntoIterator<Item=T>) -> Box<[T]> { use std::iter::FromIterator; Box::<[T]>::from_iter(iter) }
  let molar_masses = from_iter(species.values().map(|Specie{composition,..}| composition.iter().map(|(element, &count)| real(count as f32)*standard_atomic_weights[element]).sum()));
	fn index(keys: &[&str], map: &Map<&str, u8>) -> (Vec<usize>, Vec<u8>) { (map.keys().map(|&key| keys.iter().position(|&k| k==key).expect(key)).collect(), map.values().copied().collect()) }
	let index = |map| index(&from_iter(species.keys().copied()), map);
	let reactions = from_iter(Vec::from(reactions).into_iter().map(|Reaction{equation, model}| (equation, {use Model::*; match model {
			Elementary{rate_constants} => ModelParameters::Elementary{rate_constants},
			ThreeBody{rate_constants, efficiencies} => ModelParameters::ThreeBody{rate_constants, efficiencies: species.keys().map(|specie| *efficiencies.get(specie).unwrap_or(&real(1.))).collect()},
			Falloff{k0, k_inf, efficiencies, troe} => ModelParameters::Falloff{k0, k_inf, efficiencies: species.keys().map(|specie| *efficiencies.get(specie).unwrap_or(&real(1.))).collect(), troe},
	}})));
	let reactions = from_iter(species.keys().map(|specie| from_iter(
		reactions.iter().map(|(equation, model)| (- (*equation.0.get(specie).unwrap_or(&0) as i8) + *equation.1.get(specie).unwrap_or(&0) as i8, (index(&equation.0), index(&equation.1)), model))
		.filter(|&(coefficient,_,_)| coefficient!=0) )));
	let Phase::IdealGas{state, ..} = &phases[0];
	let InitialState{pressure, volume, temperature, mole_proportions} = state;
	let mole_proportions = from_iter(species.keys().map(|specie| *mole_proportions.get(specie).unwrap_or(&real(0.))));
	let total_amount = pressure * volume / (ideal_gas_constant * temperature);
	let mole_fractions = scale(real::recip(mole_proportions.iter().sum::<real>()), mole_proportions.iter().copied());
	let mass_proportions = from_iter(mul(mole_fractions, molar_masses.iter().copied()));
	let average_molar_mass = mass_proportions.iter().sum::<real>();
	let total_mass = total_amount * average_molar_mass;

	#[derive(Clone)] struct State {
		time: real,
		temperature: real,
		mass_fractions: Box<[real]>
	}
	impl From<State> for (real, Box<[Box<[real]>]>) { fn from(s: State) -> Self { (s.time, box [box [s.temperature] as Box<[_]>, s.mass_fractions] as Box<[_]>) } }
	let mut state = State{time: real(0.), temperature: *temperature, mass_fractions: from_iter(scale(real::recip(average_molar_mass), mass_proportions.iter().copied()))};
	ui::app::run(Plot{
		keys: &[&["T"] as &[_], &from_iter(species.keys().copied())],
		values: vec!(state.clone().into()),
		source: std::iter::from_fn(move || {
			let density = total_mass / (total_amount * ideal_gas_constant) * pressure / state.temperature;
			let concentrations = from_iter(scale(density, mul(recip(&molar_masses), state.mass_fractions.iter().copied())));
			let mass_fraction_rates = from_iter(scale(real::recip(density), mul(molar_masses.iter().copied(), reactions.iter().map(|reactions| reactions.iter().map(
				|(coefficient, reactants, model)| real(*coefficient as f32)*(
					rate((&reactants.0.0,&reactants.0.1), model, state.temperature, &concentrations) -
					rate((&reactants.1.0,&reactants.1.1), model, state.temperature, &concentrations)
				)
			).sum()))));
			acc(&mut state.mass_fractions, scale(time_step, mass_fraction_rates.iter().copied()));
			let specific_heat_capacities = from_iter(species.values().map(|Specie{specific_heat,..}| specific_heat.at(state.temperature)));
			let specific_heat_capacity = dot(&state.mass_fractions, &specific_heat_capacities);
			state.temperature += time_step * dot(&specific_heat_capacities, &mass_fraction_rates) / specific_heat_capacity * temperature;
			state.time += time_step;
			Some(state.clone().into())
		})
	})?
}
