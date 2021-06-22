#![allow(non_snake_case)]
use serde::{Serialize, Deserialize};
pub use {std::boxed::Box, linear_map::LinearMap as Map, strum_macros::EnumString};

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord, EnumString, Clone, Copy)] pub enum Element { H, C, N, O, Ar }

#[derive(Serialize, Deserialize, Debug, PartialEq)] pub struct State<'t> {
	pub temperature: f64,
	pub pressure: f64,
	pub volume: f64,
	#[serde(borrow)] pub amount_proportions: Map<&'t str, f64>
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Clone)] pub struct NASA7 {
	pub temperature_ranges: Box<[f64]>,
	pub pieces: Box<[[f64; 7]]>,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Clone, Copy)] pub enum Geometry {
	Atom,
	Linear {#[serde(default,rename="polarizability_A3")] polarizability_Å3: f64, #[serde(default)] rotational_relaxation: f64},
	Nonlinear {#[serde(default,rename="polarizability_A3")] polarizability_Å3: f64, #[serde(default)] rotational_relaxation: f64, #[serde(default)] permanent_dipole_moment_Debye: f64},
}
#[derive(Serialize, Deserialize, Debug, PartialEq, Clone)] pub struct Transport {
	pub well_depth_K: f64,
	#[serde(rename="diameter_A")] pub diameter_Å: f64,
	pub geometry: Geometry,
}
#[derive(Serialize, Deserialize, Debug, PartialEq, Clone)] pub struct Specie {
	pub composition: Map<Element, u8>,
	pub thermodynamic: NASA7,
	pub transport: Transport
}

#[derive(Serialize, Deserialize, Debug, PartialEq)] pub struct RateConstant {
	#[serde(rename="A")] pub preexponential_factor: f64, // m^3/mol/s
	#[serde(rename="beta")] pub temperature_exponent: f64,
	#[serde(rename="Ea")] pub activation_energy: f64 // cal/mol
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq)] pub struct Troe { pub A: f64, pub T3: f64, pub T1: f64, pub T2: f64 }

#[derive(Serialize, Deserialize, Debug, PartialEq)] pub enum ReactionModel<'t> {
	Elementary,
	Irreversible,
	ThreeBody { #[serde(borrow)] efficiencies: Map<&'t str, f64> },
	PressureModification { #[serde(borrow)] efficiencies: Map<&'t str, f64>, k0: RateConstant },
	Falloff { #[serde(borrow)] efficiencies: Map<&'t str, f64>, k0: RateConstant, troe: Troe },
}

#[derive(Serialize, Deserialize, Debug, PartialEq)] pub struct Reaction<'t> {
	#[serde(borrow)] pub equation: [Map<&'t str, u8>; 2],
	pub rate_constant: RateConstant,
	pub model: ReactionModel<'t>,
}

#[derive(Serialize, Deserialize, Debug, PartialEq)] pub struct Model<'t> {
	#[serde(borrow)] pub species: Map<&'t str, Specie>,
	#[serde(borrow)] pub reactions: Box<[Reaction<'t>]>,
	#[serde(borrow)] pub state: State<'t>,
	pub time_step: f64,
}
