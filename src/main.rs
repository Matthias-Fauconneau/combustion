enum LengthUnit { cm };
enum MolarEnergyUnit { cal_per_mol };
struct Units { length: LengthUnit, activation_energy: MolarEnergyUnit };
enum Element { H, O, Ar };
struct State { temperature: f32, pressure: f32/*atm*/};
enum Phase {
	IdealGas {
		elements: Vec<Element>,
		species: Vec<String>,
		state: State,
	}
}
struct NASA7 {
	temperature_ranges: Vec<f32>,
	data: Vec<Vec<f32>>,
}
enum Transport {
	Atom { well_depth: f32, diameter: f32},
	Linear { well_depth: f32, diameter: f32, polarizability: f32, rotational_relaxation: f32},
	Nonlinear { well_depth: f32, diameter: f32, rotational_relaxation: f32},
}
struct Specie {
	composition: BTreeMap<Element, u8>,
	thermodynamic: NASA7,
	transport: Transport
}
struct RateConstant { A: f32, b: f32, Ea: f32 };
struct Troe { A: f32, T3: f32, T1: f32, T2: f32 };
enum Reaction {
	Elementary { equation: String, rate_constant: RateConstant },
	Multiple { equation: String, rate_constants: [RateConstant; 2] },
	ThreeBody { equation: String, rate_constant: RateConstant, efficiencies: BTreeMap<String, f32> },
	Falloff { equation: String, rate_constants: [RateConstant; 2], troe: Troe, efficiencies: BTreeMap<String, f32> },
}
struct Mechanism {
	units: Units,
	phases: Vec<Phase>,
	species: BTreeMap<String, Specie>
	reactions: Vec<Reaction>,
}

fn main() {
	let mechanism: Mechanism = ron::de::from_reader(File::open("H2+O2.ron")?)?;
	dbg!(mechanism);
}
