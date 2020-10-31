#![feature(in_band_lifetimes)]
#![feature(trait_alias)]
#![feature(box_syntax)]
#![feature(map_into_keys_values)]
#![allow(non_snake_case)]
#![feature(associated_type_bounds)]
#![feature(bindings_after_at)]
pub fn scale(s: f64, v: impl IntoIterator<Item=f64,IntoIter:'t>) -> impl Iterator<Item=f64>+'t { v.into_iter().map(move |v| s*v) }
pub fn recip(x: impl IntoIterator<Item=f64,IntoIter:'t>) -> impl Iterator<Item=f64>+'t { x.into_iter().map(|x| f64::recip(x)) }
//trait Captures<'t> {}
//impl<'t, T: ?Sized> Captures<'t> for T {}
//+Captures<'a>+Captures<'b>
//pub fn mul<'a:'t,'b:'t,'t>(a: impl IntoIterator<Item=f64,IntoIter:'a>, b: impl IntoIterator<Item=f64,IntoIter:'b>) -> impl Iterator<Item=f64>+'t {
pub fn mul(a: impl IntoIterator<Item=f64,IntoIter:'t>, b: impl IntoIterator<Item=f64,IntoIter:'t>) -> impl Iterator<Item=f64>+'t {
	a.into_iter().zip(b.into_iter()).map(|(a,b)| a*b)
}
pub fn dot(a: impl IntoIterator<Item=f64,IntoIter:'t>, b: impl IntoIterator<Item=f64,IntoIter:'t>) -> f64 { mul(a.into_iter(), b.into_iter()).sum() }

mod ron {
use serde::Deserialize;
pub use std::collections::BTreeMap as Map;
#[derive(Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)] pub enum Element { H, O, Ar }
#[derive(Deserialize, Debug)] pub struct InitialState<'t> { pub temperature: f64, pub pressure: f64, #[serde(borrow)] pub mole_proportions: Map<&'t str, f64>}
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
		ideal_gas_constant * (a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T)
	}
	pub fn specific_enthalpy(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T)
	}
	pub fn specific_entropy(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[6]+a[0]*f64::ln(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T)
	}
	pub fn b(&self, T: f64) -> f64 { // S/R - H/RT
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		a[6] - a[0] + (a[0]-1.)*f64::ln(T) + a[1]/2.*T + a[2]/6.*T*T + a[3]/12.*T*T*T + a[4]/20.*T*T*T*T - a[5]/T
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

struct Coefficient { forward_rate_coefficient: f64, efficiency: f64 }
impl Model {
fn forward_rate_coefficient(&self, T: f64, concentrations: &[f64]) -> Coefficient {
	match self {
		Self::Elementary{rate_constants} => Coefficient{forward_rate_coefficient: rate_constants.iter().map(|rate_constant| arrhenius(rate_constant, T)).sum::<f64>(), efficiency: 1.},
		Self::ThreeBody{rate_constants, efficiencies} => {
			let third_body = dot(efficiencies.iter().copied(), concentrations.iter().copied());
			assert!(third_body> 0., "{} {:?} {:?}", third_body, efficiencies, concentrations);
			Coefficient{forward_rate_coefficient: rate_constants.iter().map(|rate_constant| arrhenius(rate_constant, T)).sum::<f64>(), efficiency: third_body}
		},
		Self::Falloff{k_inf, k0, efficiencies, troe: Troe{A, T3, T1, T2}} => {
			let k_inf = arrhenius(k_inf, T);
			let Pr = arrhenius(k0, T) / k_inf * dot(efficiencies.iter().copied(), concentrations.iter().copied());
			let Fcent = (1.-A)*f64::exp(-T/T3)+A*f64::exp(-T/T1)+f64::exp(-T2/T);
			let log10Fcent = f64::log10(Fcent);
			let C = -0.4-0.67*log10Fcent;
			/*let N = 0.75-1.27*log10Fcent;
			let f1 = (f64::log10(Pr) + C)/(N-0.14*(f64::log10(Pr) + C));
			let F = num::exp10(log10Fcent/(1.+f1*f1));*/
			let ATroe = f64::log10(Pr) + C;
			let BTroe = 0.806 - 1.1762*log10Fcent - 0.14*f64::log10(Pr);
			let F = f64::powf(Fcent, f64::recip(1.+num::sq(ATroe/BTroe)));
			let pressure_modification = Pr / (1.+Pr) * F; // Chemically activated bimolecular reaction
			assert!(pressure_modification > 0.);
			Coefficient{forward_rate_coefficient: k_inf, efficiency: pressure_modification}
		}
	}
}
}

pub struct Reaction {
	pub equation: [(Box<[usize]>, Box<[u8]>); 2],
	pub model: Model,
	specie_net_coefficients: Box<[f64]>
}

pub struct Rate { pub forward_base_rate: f64, pub reverse_base_rate: f64, pub efficiency: f64 }
impl Reaction {
	pub fn rate(&self, temperature: f64, B: &[f64], concentrations: &[f64]) -> Rate {
		let Reaction{equation, model, specie_net_coefficients} = self;
		let atmospheric_pressure = 101325.;
		let equilibrium_constant = f64::powf(atmospheric_pressure / ideal_gas_constant, specie_net_coefficients.iter().sum()) * f64::exp(specie_net_coefficients.iter().zip(B.iter()).map(|(net_coefficient, b)| net_coefficient * b).sum());
		let Coefficient{forward_rate_coefficient, efficiency} = model.forward_rate_coefficient(temperature, &concentrations);
		let reverse_rate_coefficient = forward_rate_coefficient / equilibrium_constant;
		let [forward_base_rate, reverse_base_rate] = iter::array::Iterator::collect::<[_;2]>(equation.iter().zip(&[forward_rate_coefficient, reverse_rate_coefficient]).map(|((species, coefficients), rate)|
			rate * species.iter().zip(coefficients.iter()).map(|(&specie, &coefficient)| concentrations[specie].powi(coefficient as i32)).product::<f64>()
		));
		//dbg!(equilibrium_constant, forward_rate_coefficient, efficiency, reverse_rate_coefficient, forward_base_rate, reverse_base_rate);
		Rate{forward_base_rate, reverse_base_rate, efficiency}
	}
}

pub struct System {
	pub molar_masses: Box<[f64]>,
	pub thermodynamics: Box<[NASA7]>,
	pub reactions: Box<[Reaction]>,
	pub average_molar_mass: f64,
	time_step: f64,
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
	let InitialState{mole_proportions, temperature, pressure, ..} = state;
	let mole_proportions = from_iter(species.iter().map(|specie| *mole_proportions.get(specie).unwrap_or(&0.)));
	let mole_fractions = from_iter(scale(f64::recip(mole_proportions.iter().sum::<f64>()), mole_proportions.iter().copied()));
	let mass_proportions = from_iter(mul(mole_fractions.iter().copied(), molar_masses.iter().copied()));
	let mass_fractions = from_iter(scale(f64::recip(mass_proportions.iter().sum::<f64>()), mass_proportions.iter().copied()));
	let average_molar_mass = f64::recip(mul(mass_fractions.iter().copied(), recip(molar_masses.iter().copied())).sum());
	{
		let density = average_molar_mass * pressure / (ideal_gas_constant * temperature);
		dbg!(density);
		let specific_enthalpy = average_molar_mass * mul(mass_fractions.iter().copied(), thermodynamics.iter().zip(molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_enthalpy(temperature) / molar_mass)).sum::<f64>();
		dbg!(specific_enthalpy);
		let specific_entropy = average_molar_mass * mul(mass_fractions.iter().copied(), thermodynamics.iter().zip(molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_entropy(temperature) / molar_mass)).sum::<f64>();
		dbg!(specific_entropy);
		let specific_heat_capacity = average_molar_mass * mul(mass_fractions.iter().copied(), thermodynamics.iter().zip(molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_heat_capacity(temperature) / molar_mass)).sum::<f64>();
		dbg!(specific_heat_capacity);
	}
	Self{
		species,
		system: System{molar_masses, thermodynamics, reactions, average_molar_mass, time_step, pressure},
		state: State{time: 0., temperature, mass_fractions}
	}
}
}

impl State {
	pub fn step(&mut self, system: &System) {
		let density = system.average_molar_mass * system.pressure / (ideal_gas_constant * self.temperature);
		let B = from_iter(system.thermodynamics.iter().map(|thermodynamic| thermodynamic.b(self.temperature)));
		let concentrations = from_iter(scale(density, mul(recip(system.molar_masses.iter().copied()), self.mass_fractions.iter().copied())));
		let specie_count = concentrations.len();
		let mut production_rates = vec!(0.; specie_count-1); // Skips most abundant specie (last index) (will be deduced from mass conservation)
		for reaction@Reaction{specie_net_coefficients,..} in system.reactions.iter() {
			let Rate{forward_base_rate, reverse_base_rate, efficiency} = reaction.rate(self.temperature, &B, &concentrations);
			let base_rate = forward_base_rate - reverse_base_rate;
			let net_rate = efficiency * base_rate;
			for (specie, net_coefficient) in specie_net_coefficients[..specie_count-1].iter().enumerate() { production_rates[specie] += net_coefficient * net_rate; }
		}
		let mass_fraction_rates = from_iter(scale(f64::recip(density), mul(system.molar_masses.iter().copied(), production_rates.iter().copied())));
		for (state, rate) in self.mass_fractions[..specie_count-1].iter_mut().zip(mass_fraction_rates.iter()) { *state = 0f64.max(*state + system.time_step * rate); }
		self.mass_fractions[specie_count-1] = 0f64.max(1. - self.mass_fractions[..specie_count-1].iter().sum::<f64>()); // Enforces mass conservation constraint by rescaling most abundant specie (last index)
		for &mass_fraction in self.mass_fractions.iter() { assert!(mass_fraction >= 0. && mass_fraction < 1.,"{} {:?}", mass_fraction, &self.mass_fractions); }
		let specific_heat_capacities = from_iter(system.thermodynamics.iter().zip(system.molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_heat_capacity(self.temperature) / molar_mass));
		let specific_heat_capacity = dot(self.mass_fractions.iter().copied(), specific_heat_capacities.iter().copied());
		let enthalpies = from_iter(system.thermodynamics.iter().zip(system.molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_enthalpy(self.temperature) / molar_mass));
		self.temperature -= system.time_step * f64::recip(specific_heat_capacity) * dot(enthalpies.iter().copied(), mass_fraction_rates.iter().copied()) ;
		self.time += system.time_step;
	}
}

impl From<State> for (f64, Box<[Box<[f64]>]>) { fn from(s: State) -> Self { (s.time*1e9/*ns*/, box [box [s.temperature] as Box<[_]>, s.mass_fractions] as Box<[_]>) } }
