#![allow(incomplete_features)]#![feature(const_generics, const_evaluatable_checked, type_ascription, non_ascii_idents,in_band_lifetimes,once_cell,array_map,map_into_keys_values,bindings_after_at,destructuring_assignment)]
#![allow(non_snake_case,confusable_idents,mixed_script_confusables,non_upper_case_globals)]
pub mod ron;
use {iter::{array_from_iter as from_iter, into::{Enumerate, IntoChain, IntoMap}, eval, vec::{eval, Dot, Prefix, Suffix}}, self::ron::{Map, Element, Troe}};

#[derive(Debug)] pub struct NASA7(pub [[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure : f64 = 101325.; // 1 atm
	const T_split : f64 = 1000.;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn dimensionless_specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T }
	pub fn dimensionless_specific_enthalpy_T(&self, T: f64) -> f64 { let a = self.a(T); a[5]/T+a[0]+a[1]/2.*T+a[2]/3.*T*T+a[3]/4.*T*T*T+a[4]/5.*T*T*T*T }
	pub fn dimensionless_specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); a[6]+a[0]*f64::ln(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T }
}

pub const ideal_gas_constant : f64 = 8.31446261815324; // J⋅K−1⋅mol−1

#[derive(Debug, Clone, Copy)] pub struct RateConstant {
	pub log_preexponential_factor: f64,
	pub temperature_exponent: f64,
	pub activation_temperature: f64
}

impl From<ron::RateConstant> for RateConstant { fn from(ron::RateConstant{preexponential_factor, temperature_exponent, activation_energy}: ron::RateConstant) -> Self {
	Self{log_preexponential_factor: f64::ln(preexponential_factor), temperature_exponent, activation_temperature: activation_energy*4.184/ideal_gas_constant}
}}

pub fn log_arrhenius(&RateConstant{log_preexponential_factor, temperature_exponent, activation_temperature}: &RateConstant, T: f64) -> f64 {
	log_preexponential_factor + temperature_exponent*f64::ln(T) - activation_temperature*(1./T)
}

#[derive(Debug, Clone, Copy)] pub enum Model<const S: usize> {
	Elementary,
	ThreeBody { efficiencies: [f64; S] },
	Falloff { efficiencies: [f64; S], k0: RateConstant, troe: Troe },
}

impl<const S: usize> Model<S> {
pub fn efficiency(&self, T: f64, concentrations: &[f64; S], log_k_inf: f64) -> f64 {
	match self {
		Self::Elementary => 1.,
		Self::ThreeBody{efficiencies} => efficiencies.dot(concentrations),
		Self::Falloff{efficiencies, k0, troe: Troe{A, T3, T1, T2}} => {
			let Pr = efficiencies.dot(concentrations) * f64::exp(log_arrhenius(k0, T) - log_k_inf); // [k0/kinf] = [C] (m3/mol)
			let Fcent = (1.-A)*f64::exp(-T/T3)+A*f64::exp(-T/T1)+f64::exp(-T2/T);
			let log10Fcent = f64::log10(Fcent);
			let C = -0.4-0.67*log10Fcent;
			let N = 0.75-1.27*log10Fcent;
			let log10PrC = f64::log10(Pr) + C;
			let f1 = log10PrC/(N-0.14*log10PrC);
			let F = num::exp10(log10Fcent/(1.+f1*f1));
			Pr / (1.+Pr) * F
		}
	}
}
}

#[derive(Clone, Copy)] pub struct Reaction<const S: usize> where [(); S-1]: {
	pub reactants: [f64; S-1],
	//Σνf: f64,
	pub products: [f64; S-1],
	//Σνr: f64,
	pub net: [f64; S-1],
	pub Σnet: f64,
	pub rate_constant: RateConstant,
	pub model: Model<S>,
}

pub struct System<const S: usize> where [(); S-1]: {
	pub molar_masses: [f64; S],
	pub thermodynamics: [NASA7; S],
	pub reactions: Box<[Reaction<S>]>,
}

impl<const S: usize> System<S> where [(); S-1]:, [(); 2+S-1]: {
	pub fn state(P: f64, y: &[f64; 2+S-1]) -> (f64, f64, [f64; S]) {
		let (T, V, n) = (y[0], y[1], y.suffix());
		let concentrations : [_; S-1] = eval(n, |n| (n / V)/*.max(0.)*/); // Skips most abundant specie (last index) (will be deduced from conservation)
		let C = P / (ideal_gas_constant * T);
		let Ca = C - concentrations.iter().sum::<f64>();
		let concentrations = from_iter(concentrations.chain([Ca]));
		(T, V, concentrations)
	}
	pub fn dt(&self, P: f64, y: &[f64; 2+S-1]) -> [f64; 2+S-1] {
		let Self{thermodynamics: species, reactions, molar_masses: W, ..} = self;
		let (T, V, concentrations) = Self::state(P, y);
		let logP0_RT = f64::ln(NASA7::reference_pressure/ideal_gas_constant) - f64::ln(T);
		let ref H_T = eval(species, |s| s.dimensionless_specific_enthalpy_T(T));
		let ref G = eval!(species, H_T; |s, h_T| h_T - s.dimensionless_specific_entropy(T)); // (H-TS)/RT
		let log_concentrations = eval(concentrations, f64::ln);
		let ref mut dtω = [0.; S-1];
		for Reaction{reactants, products, rate_constant, model, net, Σnet, ..} in reactions.iter() {
			let log_kf = log_arrhenius(rate_constant, T);
			let Rf = f64::exp(reactants.dot(log_concentrations) + log_kf);
			let log_equilibrium_constant = -net.dot(G) + Σnet*logP0_RT;
			let Rr = f64::exp(products.dot(log_concentrations) + log_kf - log_equilibrium_constant);
			let net_rate = model.efficiency(T, &concentrations, log_kf) * (Rf - Rr);
			for (specie, ν) in net.enumerate() { dtω[specie] += ν * net_rate; }
		}
		let ref dtω = *dtω;

		let Cp = species.iter().map(|s| s.dimensionless_specific_heat_capacity(T));
		let rcp_ΣCCp = 1./concentrations.dot(Cp);
		let dtT_T = - rcp_ΣCCp * dtω.dot(H_T); // R/RT
		let dtE = W.map(|w| 1.-w/W[S-1]).dot(dtω);
		let dtV = V * (dtT_T + T * ideal_gas_constant / P * dtE);
		let dtn = eval(dtω, |dtω| V*dtω);
		from_iter([dtT_T*T, dtV].chain(dtn))
	}
}

#[derive(Clone)] pub struct State<const S: usize> {
	pub temperature: f64,
	//pub volume: f64,
	pub amounts: [f64; S]
}

pub struct Simulation<'t, const S: usize> where [(); S-1]: {
	pub species: [&'t str; S],
	pub system: System<S>,
	pub time_step: f64,
	pub pressure: f64,
	pub volume: f64,
	pub state: State<S>
}

use std::lazy::SyncLazy;
pub static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	::ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}").unwrap().into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

impl<const S: usize> Simulation<'t, S> where [(); S-1]: {
	pub const species_len : usize = S;
	//pub const state_vector_len : usize = 2/*T, V*/+Self::species_len-1/*Skips most abundant specie (last index) (will be deduced from conservation)*/;
	//pub type Vec = [f64; 2/*T, V*/+S-1/*Skips most abundant specie (last index) (will be deduced from conservation)*/];

	pub fn new(system: &'b [u8]) -> ::ron::Result<Self> where 'b: 't, [(); S]: {
		let ron::System{species: species_data, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;
		let ron::Phase::IdealGas{species, state, ..} = phases.into_vec().into_iter().next().unwrap();
		use std::convert::TryInto;
		let species : [_; S] = species.as_ref().try_into().unwrap();
		let molar_masses = eval(species, |s| species_data[s].composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
		let thermodynamics = eval(species, |s| { let ron::Specie{thermodynamic: ron::NASA7{temperature_ranges, pieces},..} = &species_data[s]; match temperature_ranges[..] {
			[_,Tsplit,_] if Tsplit == NASA7::T_split => NASA7(pieces[..].try_into().unwrap()),
			[min, max] if min < NASA7::T_split && NASA7::T_split < max => NASA7([pieces[0]; 2]),
			ref ranges => panic!("{:?}", ranges),
		}});
		#[allow(unused_variables)]
		let reactions = iter::into::Collect::collect(reactions.map(|self::ron::Reaction{ref equation, rate_constant, model}| {
			let [reactants, products] = eval(equation, |e| eval(species.prefix(), |s| *e.get(s).unwrap_or(&0) as f64));
			//let [reactants, products] = eval(equation, |e| eval(species[..S-1].try_into().unwrap():[_;S-1], |s| *e.get(s).unwrap_or(&0) as f64));
			//let [reactants, products] = [[0.; S-1]; 2];
			let net = [0.; S-1]; //iter::vec::Sub::sub(&products, &reactants);
			let Σnet = 0.; //net.iter().sum();
			Reaction{
				reactants, products, net, Σnet,
				rate_constant: rate_constant.into(),
				model: {use self::ron::Model::*; match model {
					Elementary => Model::Elementary,
					ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: eval(&species, |specie| *efficiencies.get(specie).unwrap_or(&1.))},
					Falloff{efficiencies, k0, troe} => Model::Falloff{efficiencies: eval(&species, |specie| *efficiencies.get(specie).unwrap_or(&1.)), k0: k0.into(), troe},
				}},
			}
		}));

		let ron::InitialState{temperature, pressure, mole_proportions, volume} = state;
		let mole_proportions = eval(&species, |specie| *mole_proportions.get(specie).unwrap_or(&0.));
		let amount = pressure * volume / (ideal_gas_constant * temperature);
		let amounts = eval(&mole_proportions, |mole_proportion| amount/mole_proportions.iter().sum::<f64>() * mole_proportion);

		Ok(Self{
			species,
			system: System{molar_masses, thermodynamics, reactions},
			time_step, pressure, volume,
			state: State{temperature, /*volume,*/ amounts}
		})
	}
}
