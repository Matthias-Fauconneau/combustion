#![feature(in_band_lifetimes)]
#![feature(trait_alias)]
#![feature(box_syntax)]
#![feature(map_into_keys_values)]
#![allow(non_snake_case)]

pub fn scale(s: f64, v: impl Iterator<Item=f64>+'t) -> impl Iterator<Item=f64>+'t { v.map(move |v| s*v) }
pub fn recip(x: &'t [f64]) -> impl Iterator<Item=f64>+'t { x.iter().map(|&x| f64::recip(x)) }
pub fn mul(a: impl Iterator<Item=f64>+'t, b: impl Iterator<Item=f64>+'t) -> impl Iterator<Item=f64>+'t { a.zip(b).map(|(a,b)| a*b) }
pub fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b.iter()).map(|(&a,&b)| a*b).sum() }
pub fn acc(a: &mut [f64], b: impl Iterator<Item=f64>) { for (a,b) in a.iter_mut().zip(b) { *a += b; } }

mod ron {
use serde::Deserialize;
pub use std::collections::BTreeMap as Map;
#[derive(Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)] pub enum Element { H, O, Ar }
#[derive(Deserialize, Debug)] pub struct InitialState<'t> { pub temperature: f64, pub pressure: f64, pub volume: f64, #[serde(borrow)] pub mole_proportions: Map<&'t str, f64>}
#[derive(Deserialize, Debug)] pub enum Phase<'t> {
	IdealGas {
		elements: Box<[Element]>,
		species: Box<[&'t str]>,
		#[serde(borrow)] state: InitialState<'t>,
	}
}
#[derive(Deserialize, Debug)] pub struct NASA7 {
	pub temperature_ranges: Box<[f64]>,
	pub coefficients: Box<[Box<[f64]>]>,
}

#[derive(Deserialize, Debug)] enum Transport {
	Atom { well_depth: f64, diameter: f64},
	Linear { well_depth: f64, diameter: f64, polarizability: f64, rotational_relaxation: f64},
	Nonlinear { well_depth: f64, diameter: f64, rotational_relaxation: f64},
}
#[derive(Deserialize, Debug)] pub struct Specie {
	pub composition: Map<Element, u8>,
	pub thermodynamic: NASA7,
	transport: Transport
}

#[derive(Deserialize, Debug)] pub struct RateConstant {
	#[serde(rename="A")] pub preexponential_factor: f64,
	#[serde(rename="b")] pub temperature_exponent: f64,
	#[serde(rename="Ea")] pub activation_energy: f64
}

#[derive(Deserialize, Debug)] pub struct Troe { pub A: f64, pub T3: f64, pub T1: f64, pub T2: f64 }

#[derive(Deserialize, Debug)] pub enum Model<'t> {
	Elementary { rate_constants: Box<[RateConstant]> },
	ThreeBody { rate_constants: Box<[RateConstant]>, #[serde(borrow)] efficiencies: Map<&'t str, f64> },
	Falloff { k0: RateConstant, k_inf: RateConstant, #[serde(borrow)] efficiencies: Map<&'t str, f64>, troe: Troe },
}

#[derive(Deserialize, Debug)] pub struct Reaction<'t> {
	#[serde(borrow)] pub equation: [Map<&'t str, u8>; 2],
	pub model: Model<'t>,
}

#[derive(Deserialize, Debug)] pub struct System<'t> {
	pub time_step: f64,
	#[serde(borrow)] pub phases: Box<[Phase<'t>]>,
	#[serde(borrow)] pub species: Map<&'t str, Specie>,
	#[serde(borrow)] pub reactions: Box<[Reaction<'t>]>,
}
}

#[allow(non_upper_case_globals)] pub const ideal_gas_constant : f64 = 8.31446261815324;

use self::ron::*;

impl NASA7 {
	pub fn specific_heat_capacity(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		let T = T as f64;
		(ideal_gas_constant * (a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T)) as f64
	}
	pub fn specific_enthalpy(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		let T = T as f64;
		(ideal_gas_constant * (a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T)) as f64
	}
	pub fn specific_entropy(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		let T = T as f64;
		(ideal_gas_constant * (a[6]+a[0]*f64::ln(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T)) as f64
	}
}

pub fn arrhenius(&RateConstant{preexponential_factor, temperature_exponent, activation_energy}: &RateConstant, temperature: f64) -> f64 {
	preexponential_factor*temperature.powf(temperature_exponent)*f64::exp(-activation_energy/(ideal_gas_constant/4.184*temperature))
}

#[derive(Debug)] pub enum Model {
	Elementary { rate_constants: Box<[RateConstant]> },
	ThreeBody { rate_constants: Box<[RateConstant]>, efficiencies: Box<[f64]> },
	Falloff { k0: RateConstant, k_inf: RateConstant, efficiencies: Box<[f64]>, troe: Troe },
}

pub fn rate(reactants: (&[usize], &[u8]), model: &Model, T: f64, concentrations: &[f64]) -> f64 {
	let reactants : f64 = reactants.0.iter().zip(reactants.1.iter()).map(|(&specie, &coefficient)| concentrations[specie].powi(coefficient as i32)).product();
	use self::Model::*;
	reactants * match model {
		Elementary{rate_constants} => rate_constants.iter().map(|rate_constant| arrhenius(rate_constant, T)).sum::<f64>(),
		ThreeBody{rate_constants, efficiencies} => rate_constants.iter().map(|rate_constant| arrhenius(rate_constant, T) * reactants).sum::<f64>() * dot(efficiencies, concentrations),
		Falloff{k_inf, k0, efficiencies, troe: Troe{A, T3, T1, T2}} => {
			let k_inf = arrhenius(k_inf, T);
			let Pr = arrhenius(k0, T) / k_inf * dot(efficiencies, concentrations);
			let Fcent = (1.-A)*f64::exp(-T/T3)+A*f64::exp(-T/T1)+f64::exp(-T2/T);
			let log10Fcent = f64::log10(Fcent);
			let C = -0.4-0.67*log10Fcent;
			let N = 0.75-1.27*log10Fcent;
			let f1 = (f64::log10(Pr) + C)/(N-0.14*(f64::log10(Pr) + C));
			let F = num::exp10(log10Fcent/(1.+f1*f1));
			k_inf * Pr / (1.+Pr) * F
		}
	}
}

pub struct Reaction {
	pub equation: [(Box<[usize]>, Box<[u8]>); 2],
	pub model: Model,
	specie_net_coefficients: Box<[f64]>
}

pub struct System {
	pub molar_masses: Box<[f64]>,
	pub thermodynamics: Box<[NASA7]>,
	pub reactions: Box<[Reaction]>,
	time_step: f64,
	pub amount: f64,
	pub mass: f64,
	pub pressure: f64,
}

#[derive(Clone)] pub struct State {
	pub time: f64,
	pub temperature: f64,
	pub mass_fractions: Box<[f64]>
}

pub struct Simulation<'t> {
	pub species: Box<[&'t str]>,
	pub system: System,
	pub state: State
}

use iter::from_iter;

impl Simulation<'t> {
#[fehler::throws(anyhow::Error)] pub fn new(system: &'b [u8]) -> Self where 'b: 't {
	let ron::System{species, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;

	let standard_atomic_weights : Map<Element, f64> = ::ron::de::from_str("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}")?;
	let standard_atomic_weights : Map<_, f64> = standard_atomic_weights.into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect();

	let specie_names = from_iter(species.keys().copied());
	let molar_masses = from_iter(species.values().map(|Specie{composition,..}| composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum()));
	let thermodynamics = from_iter(species.into_values().map(|Specie{thermodynamic,..}| thermodynamic));
	let species = specie_names;

	let reactions = from_iter(Vec::from(reactions).into_iter().map(|self::ron::Reaction{equation, model}| Reaction{
		equation: iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|e| (e.keys().map(|&key| species.iter().position(|&k| k==key).expect(key)).collect(), e.values().copied().collect()))),
		model: {use self::ron::Model::*; match model {
			Elementary{rate_constants} => Model::Elementary{rate_constants},
			ThreeBody{rate_constants, efficiencies} => Model::ThreeBody{rate_constants, efficiencies: species.iter().map(|specie| *efficiencies.get(specie).unwrap_or(&1.)).collect()},
			Falloff{k0, k_inf, efficiencies, troe} => Model::Falloff{k0, k_inf, efficiencies: species.iter().map(|specie| *efficiencies.get(specie).unwrap_or(&1.)).collect(), troe},
		}},
		specie_net_coefficients: from_iter(species.iter().map(|specie| {let [left, right] = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|e| *e.get(specie).unwrap_or(&0) as i8)); (- left + right) as f64})),
	}));

	let Phase::IdealGas{state, ..} = Vec::from(phases).into_iter().next().unwrap();
	let InitialState{mole_proportions, temperature, pressure, volume} = state;
	let mole_proportions = from_iter(species.iter().map(|specie| *mole_proportions.get(specie).unwrap_or(&0.)));
	let mole_fractions = from_iter(scale(f64::recip(mole_proportions.iter().sum::<f64>()), mole_proportions.iter().copied()));
	let mass_proportions = from_iter(mul(mole_fractions.iter().copied(), molar_masses.iter().copied()));
	let average_molar_mass = mass_proportions.iter().sum::<f64>();
	let mass_fractions = from_iter(scale(f64::recip(average_molar_mass), mass_proportions.iter().copied()));
	let specific_enthalpy : f64 = mul(mole_fractions.iter().copied(), thermodynamics.iter().map(|thermodynamic| thermodynamic.specific_enthalpy(temperature))).sum();
	dbg!(specific_enthalpy);
	let specific_entropy : f64 = mul(mole_fractions.iter().copied(), thermodynamics.iter().map(|thermodynamic| thermodynamic.specific_entropy(temperature))).sum();
	dbg!(specific_entropy);
	let specific_heat_capacity : f64 = mul(mole_fractions.iter().copied(), thermodynamics.iter().map(|thermodynamic| thermodynamic.specific_heat_capacity(temperature))).sum();
	dbg!(specific_heat_capacity);
	let amount = pressure * volume / (ideal_gas_constant * temperature);
	let mass = amount * average_molar_mass;

	Self{
		species,
		system: System{molar_masses, thermodynamics, reactions, time_step, amount, mass, pressure},
		state: State{time: 0., temperature, mass_fractions}
	}
}
}

impl State {
	pub fn step(&mut self, system: &System) {
			let density = system.mass / (system.amount * ideal_gas_constant) * system.pressure / self.temperature;
			let concentrations = from_iter(scale(density, mul(recip(&system.molar_masses), self.mass_fractions.iter().copied())));
			let specie_count = concentrations.len();
			let mut concentration_rates = vec!(0.; specie_count-1); // Skips most abundant specie (last index) (will be deduced from mass conservation)
			for Reaction{equation, model, specie_net_coefficients} in system.reactions.iter() {
					let v = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|(species, coefficients)| rate((&species, &coefficients), model, self.temperature, &concentrations)));
					let v = v[0] - v[1];
					for (specie, net_coefficient) in specie_net_coefficients[..specie_count-1].iter().enumerate() { concentration_rates[specie] += net_coefficient * v; }
			}
			let mass_fraction_rates = from_iter(scale(f64::recip(density), mul(system.molar_masses.iter().copied(), concentration_rates.iter().copied())));
			acc(&mut self.mass_fractions, scale(system.time_step, mass_fraction_rates.iter().copied()));
			self.mass_fractions[specie_count-1] = 1. - self.mass_fractions[..specie_count-1].iter().sum::<f64>(); // Enforces mass conservation constraint by rescaling most abundant specie (last index)
			let specific_heat_capacities = from_iter(system.thermodynamics.iter().map(|thermodynamic| thermodynamic.specific_heat_capacity(self.temperature)));
			let specific_heat_capacity = dot(&self.mass_fractions, &specific_heat_capacities);
			self.temperature += system.time_step * dot(&specific_heat_capacities, &mass_fraction_rates) / specific_heat_capacity * self.temperature;
			self.time += system.time_step;
		}
}

impl From<State> for (f64, Box<[Box<[f64]>]>) { fn from(s: State) -> Self { (s.time*1e9/*ns*/, box [box [s.temperature] as Box<[_]>, s.mass_fractions] as Box<[_]>) } }
