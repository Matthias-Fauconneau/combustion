#![allow(non_camel_case_types, non_snake_case)]
use serde::Deserialize;
#[derive(Deserialize, Debug)] enum LengthUnit { cm }
#[derive(Deserialize, Debug)] enum MolarEnergyUnit { cal_per_mol }
#[derive(Deserialize, Debug)] struct Units { length: LengthUnit, activation_energy: MolarEnergyUnit }
#[derive(Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)] enum Element { H, O, Ar }
#[derive(Deserialize, Debug)] struct State { temperature: f32, pressure: f32/*atm*/}
#[derive(Deserialize, Debug)] enum Phase {
	IdealGas {
		elements: Vec<Element>,
		species: Vec<String>,
		state: State,
	}
}
#[derive(Deserialize, Debug)] struct NASA7 {
	temperature_ranges: Vec<f32>,
	data: Vec<Vec<f32>>,
}
#[derive(Deserialize, Debug)] enum Transport {
	Atom { well_depth: f32, diameter: f32},
	Linear { well_depth: f32, diameter: f32, polarizability: f32, rotational_relaxation: f32},
	Nonlinear { well_depth: f32, diameter: f32, rotational_relaxation: f32},
}
use std::collections::BTreeMap as Map;
#[derive(Deserialize, Debug)] struct Specie {
	composition: Map<Element, u8>,
	thermodynamic: NASA7,
	transport: Transport
}
#[derive(Deserialize, Debug)] struct RateConstant { A: f32, b: f32, Ea: f32 }
#[derive(Deserialize, Debug)] struct Troe { A: f32, T3: f32, T1: f32, T2: f32 }
#[derive(Deserialize, Debug)] enum Reaction {
	Elementary { equation: String, rate_constant: RateConstant },
	Multiple { equation: String, rate_constants: Vec<RateConstant> },
	ThreeBody { equation: String, rate_constant: RateConstant, efficiencies: Map<String, f32> },
	Falloff { equation: String, rate_constants: Vec<RateConstant>, troe: Troe, efficiencies: Map<String, f32> },
}
#[derive(Deserialize, Debug)] struct Mechanism {
	units: Units,
	phases: Vec<Phase>,
	species: Map<String, Specie>,
	reactions: Vec<Reaction>,
}

#[fehler::throws(anyhow::Error)] fn main() {
	let mechanism = std::fs::File::open("H2+O2.ron")?;
	let mechanism: Mechanism = ron::de::from_reader(&mechanism)?;
	dbg!(mechanism);
}
