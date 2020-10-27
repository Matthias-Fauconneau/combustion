#![feature(in_band_lifetimes)]
#![feature(trait_alias)]
#![feature(box_syntax)]
#![feature(map_into_keys_values)]
#![allow(non_snake_case)]

use num::real;
pub fn scale(s: real, v: impl Iterator<Item=real>+'t) -> impl Iterator<Item=real>+'t { v.map(move |v| s*v) }
pub fn recip(x: &'t [real]) -> impl Iterator<Item=real>+'t { x.iter().map(|&x| real::recip(x)) }
pub fn mul(a: impl Iterator<Item=real>+'t, b: impl Iterator<Item=real>+'t) -> impl Iterator<Item=real>+'t { a.zip(b).map(|(a,b)| a*b) }
pub fn dot(a: &[real], b: &[real]) -> real { a.iter().zip(b.iter()).map(|(&a,&b)| a*b).sum() }
pub fn acc(a: &mut [real], b: impl Iterator<Item=real>) { for (a,b) in a.iter_mut().zip(b) { *a += b; } }

mod ron {
use serde::Deserialize;
pub use std::collections::BTreeMap as Map;
use num::real;
#[derive(Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)] pub enum Element { H, O, Ar }
#[derive(Deserialize, Debug)] pub struct InitialState<'t> { pub temperature: real, pub pressure: real, pub volume: real, #[serde(borrow)] pub mole_proportions: Map<&'t str, real>}
#[derive(Deserialize, Debug)] pub enum Phase<'t> {
	IdealGas {
		elements: Box<[Element]>,
		species: Box<[&'t str]>,
		#[serde(borrow)] state: InitialState<'t>,
	}
}
#[derive(Deserialize, Debug)] pub struct NASA7 {
	pub temperature_ranges: Box<[real]>,
	pub coefficients: Box<[Box<[real]>]>,
}

#[derive(Deserialize, Debug)] enum Transport {
	Atom { well_depth: real, diameter: real},
	Linear { well_depth: real, diameter: real, polarizability: real, rotational_relaxation: real},
	Nonlinear { well_depth: real, diameter: real, rotational_relaxation: real},
}
#[derive(Deserialize, Debug)] pub struct Specie {
	pub composition: Map<Element, u8>,
	pub thermodynamic: NASA7,
	transport: Transport
}

#[derive(Deserialize, Debug)] pub struct RateConstant {
	#[serde(rename="A")] pub preexponential_factor: real,
	#[serde(rename="b")] pub temperature_exponent: real,
	#[serde(rename="Ea")] pub activation_energy: real
}

#[derive(Deserialize, Debug)] pub struct Troe { pub A: real, pub T3: real, pub T1: real, pub T2: real }

#[derive(Deserialize, Debug)] pub enum Model<'t> {
	Elementary { rate_constants: Box<[RateConstant]> },
	ThreeBody { rate_constants: Box<[RateConstant]>, #[serde(borrow)] efficiencies: Map<&'t str, real> },
	Falloff { k0: RateConstant, k_inf: RateConstant, #[serde(borrow)] efficiencies: Map<&'t str, real>, troe: Troe },
}

#[derive(Deserialize, Debug)] pub struct Reaction<'t> {
	#[serde(borrow)] pub equation: [Map<&'t str, u8>; 2],
	pub model: Model<'t>,
}

#[derive(Deserialize, Debug)] pub struct System<'t> {
	pub time_step: real,
	#[serde(borrow)] pub phases: Box<[Phase<'t>]>,
	#[serde(borrow)] pub species: Map<&'t str, Specie>,
	#[serde(borrow)] pub reactions: Box<[Reaction<'t>]>,
}
}

#[allow(non_upper_case_globals)] const ideal_gas_constant : real = real(8.31446261815324);

use self::ron::*;

impl NASA7 {
	fn specific_heat_capacity(&self, T: real) -> real {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T)
	}
	fn specific_enthalpy(&self, T: real) -> real {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[5]+a[0]*T+a[1]/real(2.)*T*T+a[2]/real(3.)*T*T*T+a[3]/real(4.)*T*T*T*T+a[4]/real(5.)*T*T*T*T*T)
	}
	fn specific_entropy(&self, T: real) -> real {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[6]+a[0]*real::ln(T)+a[1]*T+a[2]/real(2.)*T*T+a[3]/real(3.)*T*T*T+a[4]/real(4.)*T*T*T*T)
	}
}

pub fn arrhenius(&RateConstant{preexponential_factor, temperature_exponent, activation_energy}: &RateConstant, temperature: real) -> real {
	preexponential_factor*temperature.pow(temperature_exponent)*real::exp(-activation_energy/(ideal_gas_constant/real(4.184)*temperature))
}

#[derive(Debug)] enum Model {
	Elementary { rate_constants: Box<[RateConstant]> },
	ThreeBody { rate_constants: Box<[RateConstant]>, efficiencies: Box<[real]> },
	Falloff { k0: RateConstant, k_inf: RateConstant, efficiencies: Box<[real]>, troe: Troe },
}

fn rate(reactants: (&[usize], &[u8]), model: &Model, T: real, concentrations: &[real]) -> real {
	let reactants : real = reactants.0.iter().zip(reactants.1.iter()).map(|(&specie, &coefficient)| concentrations[specie].powi(coefficient as i32)).product();
	use self::Model::*;
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
			let F = real::exp10(log10Fcent/(real(1.)+f1*f1));
			k_inf * Pr / (real(1.)+Pr) * F
		}
	}
}

struct Reaction {
	equation: [(Box<[usize]>, Box<[u8]>); 2],
	model: Model,
	specie_net_coefficients: Box<[real]>
}

pub struct System {
	molar_masses: Box<[real]>,
	thermodynamics: Box<[NASA7]>,
	reactions: Box<[Reaction]>,
	time_step: real,
	amount: real,
	mass: real,
	pressure: real,
}

#[derive(Clone)] pub struct State {
	pub time: real,
	pub temperature: real,
	pub mass_fractions: Box<[real]>
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

	let standard_atomic_weights : Map<Element, real> = ::ron::de::from_str("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}")?;
	let standard_atomic_weights : Map<_, real> = standard_atomic_weights.into_iter().map(|(e,g)| (e, g*real(1e-3/*kg/g*/))).collect();

	let specie_names = from_iter(species.keys().copied());
	let molar_masses = from_iter(species.values().map(|Specie{composition,..}| composition.iter().map(|(element, &count)| real(count as f32)*standard_atomic_weights[element]).sum()));
	let thermodynamics = from_iter(species.into_values().map(|Specie{thermodynamic,..}| thermodynamic));
	let species = specie_names;

	let reactions = from_iter(Vec::from(reactions).into_iter().map(|self::ron::Reaction{equation, model}| Reaction{
		equation: iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|e| (e.keys().map(|&key| species.iter().position(|&k| k==key).expect(key)).collect(), e.values().copied().collect()))),
		model: {use self::ron::Model::*; match model {
			Elementary{rate_constants} => Model::Elementary{rate_constants},
			ThreeBody{rate_constants, efficiencies} => Model::ThreeBody{rate_constants, efficiencies: species.iter().map(|specie| *efficiencies.get(specie).unwrap_or(&real(1.))).collect()},
			Falloff{k0, k_inf, efficiencies, troe} => Model::Falloff{k0, k_inf, efficiencies: species.iter().map(|specie| *efficiencies.get(specie).unwrap_or(&real(1.))).collect(), troe},
		}},
		specie_net_coefficients: from_iter(species.iter().map(|specie| {let [left, right] = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|e| *e.get(specie).unwrap_or(&0) as i8)); real((- left + right) as f32)})),
	}));

	let Phase::IdealGas{state, ..} = Vec::from(phases).into_iter().next().unwrap();
	let InitialState{mole_proportions, temperature, pressure, volume} = state;
	let mole_proportions = from_iter(species.iter().map(|specie| *mole_proportions.get(specie).unwrap_or(&real(0.))));
	let mole_fractions = from_iter(scale(real::recip(mole_proportions.iter().sum::<real>()), mole_proportions.iter().copied()));
	let mass_proportions = from_iter(mul(mole_fractions.iter().copied(), molar_masses.iter().copied()));
	let average_molar_mass = mass_proportions.iter().sum::<real>();
	let mass_fractions = from_iter(scale(real::recip(average_molar_mass), mass_proportions.iter().copied()));
	let specific_enthalpy : real = mul(mole_fractions.iter().copied(), thermodynamics.iter().map(|thermodynamic| thermodynamic.specific_enthalpy(temperature))).sum();
	dbg!(specific_enthalpy);
	let specific_entropy : real = mul(mole_fractions.iter().copied(), thermodynamics.iter().map(|thermodynamic| thermodynamic.specific_entropy(temperature))).sum();
	dbg!(specific_entropy);
	let specific_heat_capacity : real = mul(mole_fractions.iter().copied(), thermodynamics.iter().map(|thermodynamic| thermodynamic.specific_heat_capacity(temperature))).sum();
	dbg!(specific_heat_capacity);
	let amount = pressure * volume / (ideal_gas_constant * temperature);
	let mass = amount * average_molar_mass;

	Self{
		species,
		system: System{molar_masses, thermodynamics, reactions, time_step, amount, mass, pressure},
		state: State{time: real(0.), temperature, mass_fractions}
	}
}
}

impl State {
	pub fn step(&mut self, system: &System) {
			let density = system.mass / (system.amount * ideal_gas_constant) * system.pressure / self.temperature;
			dbg!(density);
			let concentrations = from_iter(scale(density, mul(recip(&system.molar_masses), self.mass_fractions.iter().copied())));
			let mut concentration_rates = vec!(real(0.); concentrations.len());
			for Reaction{equation, model, specie_net_coefficients} in system.reactions.iter() {
					let v = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|e| rate((&e.0, &e.1), model, self.temperature, &concentrations)));
					let v = v[0] - v[1];
					for (specie, net_coefficient) in specie_net_coefficients.iter().enumerate() { concentration_rates[specie] += net_coefficient * v; }
			}
			let mass_fraction_rates = from_iter(scale(real::recip(density), mul(system.molar_masses.iter().copied(), concentration_rates.iter().copied())));
			acc(&mut self.mass_fractions, scale(system.time_step, mass_fraction_rates.iter().copied()));
			let specific_heat_capacities = from_iter(system.thermodynamics.iter().map(|thermodynamic| thermodynamic.specific_heat_capacity(self.temperature)));
			let specific_heat_capacity = dot(&self.mass_fractions, &specific_heat_capacities);
			self.temperature += system.time_step * dot(&specific_heat_capacities, &mass_fraction_rates) / specific_heat_capacity * self.temperature;
			self.time += system.time_step;
		}
}

impl From<State> for (real, Box<[Box<[real]>]>) { fn from(s: State) -> Self { (s.time*real(1e9)/*ns*/, box [box [s.temperature] as Box<[_]>, s.mass_fractions] as Box<[_]>) } }
