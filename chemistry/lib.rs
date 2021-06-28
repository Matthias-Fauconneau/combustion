#![feature(once_cell, in_band_lifetimes, array_map, array_methods)]
//#![feature(unboxed_closures, associated_type_bounds, once_cell, default_free_fn, fn_traits, , trait_alias)]
#![allow(non_upper_case_globals, non_snake_case)]
//#![allow(uncommon_codepoints, confusable_idents, non_upper_case_globals, non_snake_case)]

pub const kB : f64 = 1.380649e-23; // J / K
pub const NA : f64 = 6.02214076e23;
const Cm_per_Debye : f64 = 3.33564e-30; //C·m (Coulomb=A⋅s)

#[derive(PartialEq, Debug)] pub struct NASA7 {
	pub temperature_split : f64,
	pub pieces: [[f64; 7]; 2],
}

impl NASA7 {
	pub const reference_pressure : f64 = 101325. / (kB*NA);
	pub fn piece(&self, T: f64) -> &[f64; 7] { &self.pieces[if T < self.temperature_split { 0 } else { 1 }] }
	pub fn molar_heat_capacity_at_constant_pressure_R(&self, T: f64) -> f64 { let a = self.piece(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T } // /R
	pub fn enthalpy_RT(&self, T: f64) -> f64 { let a = self.piece(T); a[0]+a[1]/2.*T + a[2]/3.*T*T + a[3]/4.*T*T*T + a[4]/5.*T*T*T*T + a[5]/T } // /RT
}

use {std::lazy::SyncLazy, linear_map::LinearMap as Map, model::Element};
static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(||
	{use Element::*; [(H, 1.008), (C, 12.011), (N, 14.0067), (O, 15.999), (Ar, 39.95)]}.map(|(e,g)| (e, g/1e3/*kg/g*/)).into_iter().collect()
);

#[derive(Debug)] pub struct Species {
	pub molar_mass: Box<[f64]>,
	pub thermodynamics: Box<[NASA7]>,
	pub diameter: Box<[f64]>,
	pub well_depth_J: Box<[f64]>,
	pub polarizability: Box<[f64]>,
	pub permanent_dipole_moment: Box<[f64]>,
	pub rotational_relaxation: Box<[f64]>,
	pub internal_degrees_of_freedom: Box<[f64]>,
	pub heat_capacity_ratio: Box<[f64]>,
}

use iter::map;

impl Species {
	pub fn new(species: &Map<&'t str, model::Specie>) -> (Box<[&'t str]>, Self) {
		let molar_mass = map(species, |(_,s)| s.composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
		let thermodynamics = map(species, |(_, model::Specie{thermodynamic: model::NASA7{temperature_ranges, pieces},..})| match temperature_ranges[..] {
			[_, temperature_split, _] => NASA7{temperature_split, pieces: pieces[..].try_into().unwrap()},
			[_, _] => NASA7{temperature_split: f64::INFINITY, pieces: [pieces[0]; 2]},
			ref ranges => panic!("{:?}", ranges),
		});
		let diameter = map(species, |(_,s)| s.transport.diameter_Å*1e-10);
		let well_depth_J = map(species, |(_,s)| s.transport.well_depth_K * kB);
		use model::Geometry::*;
		let polarizability = map(species, |(_,s)| if let Linear{polarizability_Å3,..}|Nonlinear{polarizability_Å3,..} = s.transport.geometry { polarizability_Å3*1e-30 } else { 0. });
		let permanent_dipole_moment = map(species, |(_,s)|
			if let Nonlinear{permanent_dipole_moment_Debye,..} = s.transport.geometry { permanent_dipole_moment_Debye*Cm_per_Debye } else { 0. });
		let rotational_relaxation =
			map(species, |(_,s)| if let Linear{rotational_relaxation,..}|Nonlinear{rotational_relaxation,..} = s.transport.geometry { rotational_relaxation } else { 0. });
		let internal_degrees_of_freedom = map(species, |(_,s)| match s.transport.geometry { Atom => 0., Linear{..} => 1., Nonlinear{..} => 3./2. });
		let heat_capacity_ratio = map(species, |(_,s)| 1. + 2. / match s.transport.geometry { Atom => 3., Linear{..} => 5., Nonlinear{..} => 6. });
		(map(species, |(name,_)| *name), Species{molar_mass, thermodynamics, diameter, well_depth_J, polarizability, permanent_dipole_moment, rotational_relaxation, internal_degrees_of_freedom, heat_capacity_ratio})
	}
	pub fn len(&self) -> usize { self.molar_mass.len() }
}

pub struct State {
    pub temperature: f64,
    pub pressure_R: f64,
    pub volume: f64,
    pub amounts: Box<[f64]>
}

pub fn initial_state(model::Model{species, state, ..}: &model::Model<'t>) -> State {
	let model::State{temperature, pressure, volume, amount_proportions} = state;
	let species_names = map(species, |(name,_)| *name);
	for (specie,_) in amount_proportions { assert!(species_names.contains(specie)); }
	let amount_proportions = map(&*species_names, |specie| *amount_proportions.get(specie).unwrap_or(&0.));
	let pressure_R = pressure/(kB*NA);
	let temperature = *temperature; //K*: K->J
	let amount = pressure_R * volume / temperature;
	let amounts = amount_proportions.iter().map(|amount_proportion| amount * amount_proportion/amount_proportions.iter().sum::<f64>()).collect();
	State{temperature, pressure_R, volume: *volume, amounts}
}

#[derive(Clone, Copy)] pub struct RateConstant {
	pub preexponential_factor: f64,
	pub temperature_exponent: f64,
	pub activation_temperature: f64
}

impl From<&model::RateConstant> for RateConstant {
	fn from(model::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: &model::RateConstant) -> Self {
		const J_per_cal: f64 = 4.184;
		Self{preexponential_factor: *preexponential_factor, temperature_exponent: *temperature_exponent, activation_temperature: activation_energy*J_per_cal/(kB*NA)}
	}
}

pub use model::Troe;

pub enum ReactionModel {
	Elementary,
	Irreversible,
	ThreeBody { efficiencies: Box<[f64]> },
	PressureModification { efficiencies: Box<[f64]>, k0: RateConstant },
	Falloff { efficiencies: Box<[f64]>, k0: RateConstant, troe: Troe },
}

pub struct Reaction {
	pub reactants: Box<[u8]>,
	pub products: Box<[u8]>,
	pub net: Box<[i8/*; S-1*/]>,
	pub Σreactants: u8,
	pub Σproducts: u8,
	pub Σnet: i8,
	pub rate_constant: RateConstant,
	pub model: ReactionModel,
}

impl Reaction {
	pub fn new(species_names: &[&str], model::Reaction{ref equation, rate_constant, model}: &model::Reaction) -> Self {
		for side in equation { for (specie, _) in side { assert!(species_names.contains(&specie), "{}", specie) } }
		let [reactants, products] = equation.each_ref().map(|e| species_names.iter().map(|&s| *e.get(s).unwrap_or(&0)).collect::<Box<_>>());
		let net = products.into_iter().zip(reactants.into_iter()).take(species_names.len()-1).map(|(&a, &b)| a as i8 - b as i8).collect();
		let [Σreactants, Σproducts] = [reactants.iter().sum(), products.iter().sum()];
		let Σnet = Σproducts as i8 - Σreactants as i8;
		let from = |efficiencies:&Map<_,_>| map(species_names, |&specie| *efficiencies.get(specie).unwrap_or(&1.));
		Reaction{
			reactants, products, net, Σreactants, Σproducts, Σnet,
			rate_constant: rate_constant.into(),
			model: {use model::ReactionModel::*; match model {
				Elementary => ReactionModel::Elementary,
				Irreversible => ReactionModel::Irreversible,
				ThreeBody{efficiencies} => ReactionModel::ThreeBody{efficiencies: from(efficiencies)},
				PressureModification{efficiencies, k0} => ReactionModel::PressureModification{efficiencies: from(efficiencies), k0: k0.into()},
				Falloff{efficiencies, k0, troe} => ReactionModel::Falloff{efficiencies: from(efficiencies), k0: k0.into(), troe: *troe},
			}}
		}
	}
}