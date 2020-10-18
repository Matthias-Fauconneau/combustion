#![feature(in_band_lifetimes)]
#![allow(non_camel_case_types, non_snake_case,non_upper_case_globals)]

fn scale(s: f32, v: impl Iterator<Item=f32>+'t) -> impl Iterator<Item=f32>+'t { v.map(move |v| s*v) }
fn rcp(x: &'t [f32]) -> impl Iterator<Item=f32>+'t { x.iter().map(|x| 1./x) }
fn mul(a: impl Iterator<Item=f32>+'t, b: impl Iterator<Item=f32>+'t) -> impl Iterator<Item=f32>+'t { a.zip(b).map(|(a,b)| a*b) }
fn dot(a: &[f32], b: &[f32]) -> f32 { a.iter().zip(b.iter()).map(|(a,b)| a*b).sum() }
fn acc(a: &mut [f32], b: impl Iterator<Item=f32>) { for (a,b) in a.iter_mut().zip(b) { *a += b; } }

use std::collections::BTreeMap as Map;
use serde::Deserialize;
#[derive(Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)] enum Element { H, O, Ar }
#[derive(Deserialize, Debug)] struct State { temperature: f32, pressure: f32, volume: f32, mole_proportions: Map<String, f32>}
#[derive(Deserialize, Debug)] enum Phase {
	IdealGas {
		elements: Vec<Element>,
		species: Vec<String>,
		state: State,
	}
}
#[derive(Deserialize, Debug)] struct NASA7 {
	temperature_ranges: Vec<f32>,
	coefficients: Vec<Vec<f32>>,
}
impl NASA7 {
	fn at(&self, T: f32) -> f32 {
		let a = &self.coefficients[self.temperature_ranges.windows(2).position(|range| range[0] <= T && T <= range[1]).unwrap_or_else(|| panic!("{}", T))];
		ideal_gas_constant * (a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T)
	}
}

#[derive(Deserialize, Debug)] enum Transport {
	Atom { well_depth: f32, diameter: f32},
	Linear { well_depth: f32, diameter: f32, polarizability: f32, rotational_relaxation: f32},
	Nonlinear { well_depth: f32, diameter: f32, rotational_relaxation: f32},
}
#[derive(Deserialize, Debug)] struct Specie {
	composition: Map<Element, u8>,
	specific_heat: NASA7,
	transport: Transport
}

const ideal_gas_constant : f32 = 8.31446261815324;

#[derive(Deserialize, Debug)] struct RateConstant {
	#[serde(rename="A")] preexponential_factor: f32,
	#[serde(rename="b")] temperature_exponent: f32,
	#[serde(rename="Ea")] activation_energy: f32
}
fn arrhenius(&RateConstant{preexponential_factor, temperature_exponent, activation_energy}: &RateConstant, temperature: f32) -> f32 {
	preexponential_factor*f32::powf(temperature, temperature_exponent)*f32::exp(-activation_energy/(ideal_gas_constant/4.184*temperature))
}

#[derive(Deserialize, Debug)] struct Troe { A: f32, T3: f32, T1: f32, T2: f32 }
#[derive(Deserialize, Debug)] enum Model {
	Elementary { rate_constants: Vec<RateConstant> },
	ThreeBody { rate_constants: Vec<RateConstant>, efficiencies: Map<String, f32> },
	Falloff { k0: RateConstant, k_inf: RateConstant, efficiencies: Map<String, f32>, troe: Troe },
}
#[derive(Debug)] enum ModelParameters {
	Elementary { rate_constants: Vec<RateConstant> },
	ThreeBody { rate_constants: Vec<RateConstant>, efficiencies: Vec<f32> },
	Falloff { k0: RateConstant, k_inf: RateConstant, efficiencies: Vec<f32>, troe: Troe },
}

#[derive(Deserialize, Debug)] struct Reaction {
	equation: (Map<String, u8>, Map<String, u8>),
	model: Model,
}

fn rate(reactants: (&[usize], &[u8]), model: &ModelParameters, T: f32, concentrations: &[f32]) -> f32 {
	let reactants : f32 = reactants.0.iter().zip(reactants.1.iter()).map(|(&specie, &coefficient)| f32::powi(concentrations[specie], coefficient as i32)).sum();
	use ModelParameters::*;
	reactants * match model {
		Elementary{rate_constants} => rate_constants.iter().map(|rate_constant| arrhenius(rate_constant, T)).sum::<f32>(),
		ThreeBody{rate_constants, efficiencies} => rate_constants.iter().map(|rate_constant| arrhenius(rate_constant, T) * reactants).sum::<f32>() * dot(efficiencies, concentrations),
		Falloff{k_inf, k0, efficiencies, troe: Troe{A, T3, T1, T2}} => {
			let k_inf = arrhenius(k_inf, T);
			let Pr = arrhenius(k0, T) / k_inf * dot(efficiencies, concentrations);
			let Fcent = (1.-A)*f32::exp(-T/T3)+A*f32::exp(-T/T1)+f32::exp(-T2/T);
			let log10Fcent = f32::log10(Fcent);
			let C = -0.4-0.67*log10Fcent;
			let N = 0.75-1.27*log10Fcent;
			let f1 = (f32::log10(Pr) + C)/(N-0.14*(f32::log10(Pr) + C));
			let F = f32::exp(f32::ln(10.)*log10Fcent/(1.+f1*f1));
			k_inf * Pr / (1.+Pr) * F
		}
	}
}

#[derive(Deserialize, Debug)] struct Mechanism {
	time_step: f32,
	phases: Vec<Phase>,
	species: Map<String, Specie>,
	reactions: Vec<Reaction>,
}

#[fehler::throws(anyhow::Error)] fn main() {
	let standard_atomic_weights : Map<Element, f32> = ron::de::from_str("{H: 1.008, O: 15.999, Ar: 39.95}")?;
	let mechanism = std::fs::File::open("H2+O2.ron")?;
	let Mechanism{species, reactions, phases, time_step} = ron::de::from_reader(&mechanism)?;
	let molar_masses = species.values().map(|Specie{composition,..}| composition.iter().map(|(element, &count)| (count as f32)*standard_atomic_weights[element]).sum()).collect::<Vec<_>>();
	fn index(keys: &[&String], map: &Map<String, u8>) -> (Vec<usize>, Vec<u8>) { (map.keys().map(|key| keys.iter().position(|k| k==&key).expect(key)).collect(), map.values().copied().collect()) }
	let index = |map| index(&species.keys().collect::<Vec<_>>(), map);
	let reactions = reactions.into_iter().map(|Reaction{equation, model}| (equation, {use Model::*; match model {
			Elementary{rate_constants} => ModelParameters::Elementary{rate_constants},
			ThreeBody{rate_constants, efficiencies} => ModelParameters::ThreeBody{rate_constants, efficiencies: species.keys().map(|specie| *efficiencies.get(specie).unwrap_or(&1.)).collect()},
			Falloff{k0, k_inf, efficiencies, troe} => ModelParameters::Falloff{k0, k_inf, efficiencies: species.keys().map(|specie| *efficiencies.get(specie).unwrap_or(&1.)).collect(), troe},
	}})).collect::<Vec<_>>();
	let reactions = species.keys().map(|specie|
		reactions.iter().map(|(equation, model)| (- (*equation.0.get(specie).unwrap_or(&0) as i8) + *equation.1.get(specie).unwrap_or(&0) as i8, (index(&equation.0), index(&equation.1)), model))
		.filter(|&(coefficient,_,_)| coefficient!=0)
		.collect::<Vec<_>>()
	).collect::<Vec<_>>();
	let Phase::IdealGas{state, ..} = &phases[0];
	let State{pressure, volume, temperature, mole_proportions} = state;
	let mole_proportions = species.keys().map(|specie| *mole_proportions.get(specie).unwrap_or(&0.)).collect::<Vec<_>>();
	let total_amount = pressure * volume / (ideal_gas_constant * temperature);
	let mole_fractions = scale(1./mole_proportions.iter().sum::<f32>(), mole_proportions.iter().copied());
	let mass_proportions = mul(mole_fractions, molar_masses.iter().copied()).collect::<Vec<_>>();
	let average_molar_mass = mass_proportions.iter().sum::<f32>();
	let total_mass = total_amount * average_molar_mass;

	let mut mass_fractions = scale(1./average_molar_mass, mass_proportions.into_iter()).collect::<Vec<_>>();
	let mut temperature = *temperature;
	for _time in (0..100000).map(|i| i as f32*time_step) {
		let density = total_mass / (total_amount * ideal_gas_constant) * pressure / temperature;
		let concentrations = scale(density, mul(rcp(&molar_masses), mass_fractions.iter().copied())).collect::<Vec<_>>();
		let mass_fraction_rates = scale(1./density, mul(molar_masses.iter().copied(), reactions.iter().map(|reactions| reactions.iter().map(
			|(coefficient, reactants, model)| (*coefficient as f32)*(rate((&reactants.0.0,&reactants.0.1), model, temperature, &concentrations) - rate((&reactants.1.0,&reactants.1.1), model, temperature, &concentrations))
		).sum()))).collect::<Vec<_>>();
		acc(&mut mass_fractions, scale(time_step, mass_fraction_rates.iter().copied()));
		let specific_heat_capacities = species.values().map(|Specie{specific_heat,..}| specific_heat.at(temperature)).collect::<Vec<_>>();
		let specific_heat_capacity = dot(&mass_fractions, &specific_heat_capacities);
		temperature += time_step * dot(&specific_heat_capacities, &mass_fraction_rates) / specific_heat_capacity * temperature;
	}
	use itertools::Itertools; println!("{:?} {:?}", species.keys().zip(mass_fractions.iter()).format(" "), temperature);
}
