use serde::Deserialize;
pub use {std::boxed::Box, linear_map::LinearMap as Map};
#[derive(Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)] pub enum Element { H, O, C, Ar }
#[derive(Deserialize, Debug)] pub struct InitialState<'t> { pub temperature: f64, pub pressure: f64, #[serde(borrow)] pub mole_proportions: Map<&'t str, f64> }
#[derive(Deserialize, Debug)] pub enum Phase<'t> {
	IdealGas {
		elements: Box<[Element]>,
		species: Box<[&'t str]>,
		#[serde(borrow)] state: InitialState<'t>,
	}
}
#[derive(Deserialize, Debug)] pub struct NASA7 {
	pub temperature_ranges: Box<[f64]>,
	pub pieces: Box<[[f64; 7]]>,
}

#[derive(Deserialize, Debug)] pub enum Geometry {
	Atom,
	Linear {#[serde(default)] polarizability: f64, #[serde(default)] rotational_relaxation: f64},
	Nonlinear {#[serde(default)] polarizability: f64, #[serde(default)] rotational_relaxation: f64, #[serde(default)] dipole: f64},
}
#[derive(Deserialize, Debug)] pub struct Transport {
	pub well_depth: f64,
	pub diameter: f64,
	geometry: Geometry,
}
#[derive(Deserialize, Debug)] pub struct Specie {
	pub composition: Map<Element, u8>,
	pub thermodynamic: NASA7,
	pub transport: Transport
}

#[derive(Deserialize, Debug)] pub struct RateConstant {
	#[serde(rename="A")] pub preexponential_factor: f64, // m^3/mol/s
	#[serde(rename="beta")] pub temperature_exponent: f64,
	#[serde(rename="Ea")] pub activation_energy: f64 // cal/mol
}

#[derive(Deserialize, Debug, Clone, Copy)] pub struct Troe { pub A: f64, pub T3: f64, pub T1: f64, pub T2: f64 }

#[derive(Deserialize, Debug)] pub enum Model<'t> {
	Elementary,
	ThreeBody { #[serde(borrow)] efficiencies: Map<&'t str, f64> },
	PressureModification { #[serde(borrow)] efficiencies: Map<&'t str, f64>, k0: RateConstant },
	Falloff { #[serde(borrow)] efficiencies: Map<&'t str, f64>, k0: RateConstant, troe: Troe },
}

#[derive(Deserialize, Debug)] pub struct Reaction<'t> {
	#[serde(borrow)] pub equation: [Map<&'t str, u8>; 2],
	pub rate_constant: RateConstant,
	pub model: Model<'t>,
}

#[derive(Deserialize, Debug)] pub struct System<'t> {
	pub time_step: f64,
	#[serde(borrow)] pub phases: Box<[Phase<'t>]>,
	#[serde(borrow)] pub species: Map<&'t str, Specie>,
	#[serde(borrow)] pub reactions: Box<[Reaction<'t>]>,
}
