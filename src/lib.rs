#![feature(min_const_generics,non_ascii_idents,in_band_lifetimes,once_cell,array_map,map_into_keys_values,bindings_after_at,destructuring_assignment)]
#![allow(non_snake_case,confusable_idents,mixed_script_confusables,non_upper_case_globals)]
pub mod ron;
use {std::convert::TryInto,
				iter::{array_from_iter as from_iter, box_collect, into::{Collect, Enumerate, IntoChain, Zip, Find}, eval, //zip,
				vec::{eval, generate, Dot, Suffix}}};

pub use self::ron::*;

pub struct NASA7([[f64; 7]; 2]);
impl NASA7 {
	pub const reference_pressure : f64 = 101325.; // 1 atm
	const T_split : f64 = 1000.;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn dimensionless_specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T }
	pub fn dimensionless_specific_enthalpy_T(&self, T: f64) -> f64 { let a = self.a(T); a[5]/T+a[0]+a[1]/2.*T+a[2]/3.*T*T+a[3]/4.*T*T*T+a[4]/5.*T*T*T*T }
	pub fn dimensionless_specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); a[6]+a[0]*f64::ln(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T }
}

pub const ideal_gas_constant : f64 = 8.31446261815324; // J⋅K−1⋅mol−1

pub struct RateConstant {
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

pub enum Model<const S: usize> {
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

pub struct Reaction<const S: usize, const S1: usize> {
	pub equation: [Zip<Box<[usize]>, Box<[f64]>>; 2],
	pub rate_constant: RateConstant,
	pub model: Model<S>,
	pub specie_net_coefficients: [f64; S1],
	//Σνf: f64,
	//Σνr: f64,
	pub sum_net_coefficients: f64,
}

pub struct System<const S: usize, const S1: usize, const N: usize> {
	pub molar_masses: [f64; S],
	pub thermodynamics: [NASA7; S],
	pub reactions: Box<[Reaction<S,S1>]>,
}

impl<const S: usize, const S1: usize, const N: usize> System<S,S1,N> {
pub fn state(P: f64, y: &[f64; N]) -> (f64, f64, [f64; S]) {
	let (T, V, n) = (y[0], y[1], y.suffix());
	let concentrations : [_; /*S-1*/S1] = eval(n, |n| (n / V).max(0.)); // Skips most abundant specie (last index) (will be deduced from conservation) // Clamps negative values yielded by integration
	let C = P / (ideal_gas_constant * T);
	let Ca = C - concentrations.iter().sum::<f64>();
	let concentrations = from_iter(concentrations.chain([Ca]));
	(T, V, concentrations)
}
pub fn dt(&self, P: f64, y: &[f64; N]) -> [f64; N] {
	let Self{thermodynamics: species, reactions, molar_masses: W, ..} = self;
	let (T, V, concentrations) = Self::state(P, y);
	let logP0_RT = f64::ln(NASA7::reference_pressure/ideal_gas_constant) - f64::ln(T);
	let ref H_T = eval(species, |s| s.dimensionless_specific_enthalpy_T(T));
	let ref G = eval!(species, H_T; |s, h_T| h_T - s.dimensionless_specific_entropy(T)); // (H-TS)/RT
	let log_concentrations = eval(concentrations, f64::ln);
	let ref mut dtω = [0.; /*S-1*/S1];
	for Reaction{equation, rate_constant, model, specie_net_coefficients: ν, sum_net_coefficients, ..} in reactions.iter() {
		let log_equilibrium_constant = sum_net_coefficients*logP0_RT - ν.dot(G);
		let log_kf = log_arrhenius(rate_constant, T);
		let log_kr = log_kf - log_equilibrium_constant;
		use iter::into::{IntoMap, Sum};
		let [log_ΠCνf, log_ΠCνr] : [f64;2] = eval(equation, |side| side.map(|(&specie, ν)| ν*log_concentrations[specie]).sum());
		let [Rf, Rr] = [log_kf + log_ΠCνf, log_kr + log_ΠCνr].map(f64::exp);
		let net_rate = model.efficiency(T, &concentrations, log_kf) * (Rf - Rr);
		for (specie, ν) in ν.enumerate() { dtω[specie] += ν * net_rate; }
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

pub fn step(&self, _rtol: f64, _atol: f64, tmax: f64, P: f64, mut u: [f64; N]) -> [f64; N] {
	const steps : usize = 600;
	let dt = tmax/(steps as f64);
	for _ in 0..steps {
		let ref k0 = self.dt(P, &u);
		u = if true {
			eval!(&u, k0; |u, k0| u+dt*k0)
		} else {
			let a = ((),
			[0.161],
			[-0.008480655492356989,0.335480655492357],
			[2.8971530571054935, -6.359448489975075, 4.3622954328695815],
			[5.325864828439257, -11.748883564062828, 7.4955393428898365, -0.09249506636175525],
			[5.86145544294642, -12.92096931784711, 8.159367898576159, -0.071584973281401, -0.028269050394068383],
			[0.09646076681806523, 0.01, 0.4798896504144996, 1.379008574103742, -3.290069515436081, 2.324710524099774]);
			let ref k1 = self.dt(P, &eval!(&u, k0; |u, k0| u+dt*(a.1[0]*k0)));
			let ref k2 = self.dt(P, &eval!(&u, k0, k1; |u, k0, k1| u+dt*(a.2[0]*k0+a.2[1]*k1)));
			let ref k3 = self.dt(P, &eval!(&u, k0, k1, k2; |u, k0, k1, k2| u+dt*(a.3[0]*k0+a.3[1]*k1+a.3[2]*k2)));
			let ref k4 = self.dt(P, &eval!(&u, k0, k1, k2, k3; |u, k0, k1, k2, k3| u+dt*(a.4[0]*k0+a.4[1]*k1+a.4[2]*k2+a.4[3]*k3)));
			let ref k5 = self.dt(P, &eval!(&u, k0, k1, k2, k3, k4; |u, k0, k1, k2, k3, k4| u+dt*(a.5[0]*k0+a.5[1]*k1+a.5[2]*k2+a.5[3]*k3+a.5[4]*k4)));
			eval!(&u, k0, k1, k2, k3, k4, k5; |u, k0, k1, k2, k3, k4, k5| u+dt*(a.6[0]*k0+a.6[1]*k1+a.6[2]*k2+a.6[3]*k3+a.6[4]*k4+a.6[5]*k5))
		};
	}
	u
}

}

#[derive(Clone)] pub struct State<const S: usize> {
	pub temperature: f64,
	//pub volume: f64,
	pub amounts: [f64; S]
}

pub struct Simulation<'t, const S: usize, const S1: usize, const N: usize> {
	pub species: [&'t str; S],
	pub system: System<S,S1,N>,
	pub time_step: f64,
	pub pressure: f64,
	pub volume: f64,
	pub state: State<S>
}

use std::lazy::SyncLazy;
pub static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	::ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}").unwrap().into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

impl<const S: usize, const S1: usize, const N: usize> Simulation<'t, S, S1, N> {
#[fehler::throws(anyhow::Error)] pub fn new(system: &'b [u8]) -> Self where 'b: 't {
	use iter::into::{IntoMap, IntoZip};
	let ron::System{species: species_data, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;
	let Phase::IdealGas{species, state, ..} = phases.into_vec().into_iter().next().unwrap();
	let species : [_; S] = species.as_ref().try_into().unwrap();
	let molar_masses = eval(species, |s| species_data[s].composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
	let thermodynamics = eval(species, |s| { let Specie{thermodynamic: ron::NASA7{temperature_ranges, pieces},..} = &species_data[s]; match temperature_ranges[..] {
		[_,Tsplit,_] if Tsplit == NASA7::T_split => NASA7(pieces[..].try_into().unwrap()),
		[min, max] if min < NASA7::T_split && NASA7::T_split < max => NASA7([pieces[0]; 2]),
		ref ranges => panic!("{:?}", ranges),
	}});
	let reactions = reactions.map(|self::ron::Reaction{equation: ref str_equation, rate_constant, model}| {
		let equation = eval(str_equation, |e| box_collect(e.keys().map(|&key| species.iter().position(|&k| k==key).expect(key))).zip(box_collect(e.values().map(|&ν| ν as f64))));
		let specie_net_coefficients = generate(|specie| {
			let [reactant, product]:[_;2] = eval(&equation, |side| side.find(|(&s,_)| s==specie).map(|(_,&ν)| ν as i8).unwrap_or(0));
			product-reactant
		});
		let mut net = Map::new();
		for (s, ν) in specie_net_coefficients.enumerate() {
			for (element, &count) in &species_data[species[s]].composition {
				if !net.contains_key(&element) { net.insert(element, 0); }
				*net.get_mut(&element).unwrap() += ν * count as i8;
			}
		}
		let sum_net_coefficients = specie_net_coefficients.iter().sum::<i8>() as f64;
		Reaction{
			equation,
			rate_constant: rate_constant.into(),
			model: {use self::ron::Model::*; match model {
				Elementary => Model::Elementary,
				ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: eval(&species, |specie| *efficiencies.get(specie).unwrap_or(&1.))},
				Falloff{efficiencies, k0, troe} => Model::Falloff{efficiencies: eval(&species, |specie| *efficiencies.get(specie).unwrap_or(&1.)), k0: k0.into(), troe},
			}},
			specie_net_coefficients: eval(specie_net_coefficients, |ν| ν as f64),
			sum_net_coefficients,
		}
	}).collect();

	let InitialState{temperature, pressure, mole_proportions, volume} = state;
	let mole_proportions = eval(&species, |specie| *mole_proportions.get(specie).unwrap_or(&0.));
	let amount = pressure * volume / (ideal_gas_constant * temperature);
	let amounts = eval(&mole_proportions, |mole_proportion| amount/mole_proportions.iter().sum::<f64>() * mole_proportion);

	Self{
		species,
		system: System{molar_masses, thermodynamics, reactions},
		time_step, pressure, volume,
		state: State{temperature, /*volume,*/ amounts}
	}
}
}
