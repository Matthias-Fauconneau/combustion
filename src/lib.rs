#![feature(min_const_generics,non_ascii_idents,in_band_lifetimes,once_cell,array_map,map_into_keys_values,bindings_after_at,destructuring_assignment)]
#![allow(non_snake_case,confusable_idents,mixed_script_confusables,non_upper_case_globals)]
pub mod ron;
use {std::convert::TryInto, num::{norm, error},
				iter::{array_from_iter as from_iter, box_collect, into::{Collect, Enumerate, IntoChain, Zip, Find, IntoMap, IntoZip}, zip, map, eval,
				vec::{eval, generate, Dot, Suffix, Sub}}};

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
	let concentrations : [_; /*S-1*/S1] = eval(n, |n| (n / V)/*.max(0.)*/); // Skips most abundant specie (last index) (will be deduced from conservation) // Clamps negative values yielded by integration
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
	for Reaction{equation, rate_constant, model, specie_net_coefficients, sum_net_coefficients, ..} in reactions.iter() {
		let log_equilibrium_constant = sum_net_coefficients*logP0_RT - specie_net_coefficients.dot(G);
		let log_kf = log_arrhenius(rate_constant, T);
		let log_kr = log_kf - log_equilibrium_constant;
		use iter::into::Sum;
		let [log_ΠCνf, log_ΠCνr] : [f64;2] = eval(equation, |side| side.map(|(&specie, ν)| ν*log_concentrations[specie]).sum());
		let [Rf, Rr] = [log_kf + log_ΠCνf, log_kr + log_ΠCνr].map(f64::exp);
		let net_rate = model.efficiency(T, &concentrations, log_kf) * (Rf - Rr);
		for (specie, ν) in specie_net_coefficients.enumerate() { dtω[specie] += ν * net_rate; }
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

// Estimate principal eigenvector/value of dyF|y
fn power_iteration(&self, P: f64, tmax: f64, y: &[f64; N], dty: &[f64; N], v: &[f64; N]) -> ([f64; N], f64) {
	let [norm_y, norm_v] = [y,v].map(norm);
	assert!(norm_y > 0.);
	let ε = norm_y * f64::EPSILON.sqrt();
	assert!(norm_v > 0.);
	let ref mut yεv = eval!(y, v; |y, v| y + ε * v / norm_v);
	let mut ρ = 0.;
	for i in 1..=50 {
		let ref dtεv = self.dt(P, yεv).sub(dty);
		let norm_dtεv = norm(dtεv);
		assert!(norm_dtεv > 0.);
		let previous_ρ = ρ;
		ρ = norm_dtεv / ε;
		if i >= 2 && f64::abs(ρ - previous_ρ) <= 0.01*ρ.max(1./tmax) { break; }
		*yεv = eval!(y, dtεv; |y, dtεv| y + (ε / norm_dtεv) * dtεv);
	}
	(yεv.sub(y), ρ * 1.2)
}

pub fn step(&self, relative_tolerance: f64, absolute_tolerance: f64, tmax: f64, P: f64, mut u: [f64; N]) -> [f64; N] {
	let mut dtu = self.dt(P, &u);
	let (mut v, mut jacobian_spectral_radius) = self.power_iteration(P, tmax, &u, &dtu, &dtu);
	let max_stages = ((relative_tolerance / (10. * f64::EPSILON)).sqrt().round() as usize).max(2);
	let mut dt = {
		let dt = (1./jacobian_spectral_radius).min(tmax);
		let ref dtu1 = self.dt(P, &eval!(&u, &dtu; |u, dtu| u + dt * dtu));
		(dt/(dt*error(zip!(&dtu, dtu1, &u).map(|(dtu, dtu1, u):(&f64,&f64,&f64)| (dtu1 - dtu) / (absolute_tolerance + relative_tolerance * u.abs())))) / 10.).min(tmax)
	};
	let (mut previous_error, mut previous_dt) = (0., 0.);
	let mut nstep = 0;
	let mut t = 0.;
	loop {
		let stages = 1 + (1. + 1.54 * dt * jacobian_spectral_radius).sqrt().floor() as usize;
		let stages = if stages > max_stages {
			dt = (max_stages*max_stages - 1) as f64 / (1.54 * jacobian_spectral_radius);
			max_stages
		} else { stages };
		let ref u1 = {
			let w0 = 1. + 2. / (13.0 * (stages * stages) as f64);
			let sqw01 = w0*w0 - 1.;
			let arg = stages as f64 * (w0 + sqw01.sqrt()).ln();
			let w1 = arg.sinh() * sqw01 / (arg.cosh() * stages as f64 * sqw01.sqrt() - w0 * arg.sinh());
			let mut B = [1. / (4.*w0*w0); 2];
			let mu_t = w1 * B[0];
			let [ref mut u0, mut u1] = [u, eval!(&u, &dtu; |u, dtu| u + mu_t * dt * dtu)];
			let mut Z = [w0, 1.];
			let mut dZ = [1., 0.];
			let mut ddZ = [0., 0.];
			for _ in 1..stages {
				let z = 2. * w0 * Z[0] - Z[1];
				let dz = 2. * w0 * dZ[0] - dZ[1] + 2. * Z[0];
				let ddz = 2. * w0 * ddZ[0] - ddZ[1] + 4. * dZ[0];
				let b = ddz / (dz * dz);
				let gamma_t = 1. - (Z[0] * B[0]);
				let nu = - b / B[1];
				let mu = 2. * b * w0 / B[0];
				let mu_t = mu * w1 / w0;
				let ref dtu1 = self.dt(P, &u1);
				for (u0, u1, dtu1, u, dtu) in zip!(u0, &mut u1, dtu1, &u, &dtu) {
					let u0_ = *u0;
					*u0 = *u1;
					*u1 = (1.-mu-nu)*u + nu*u0_ + mu**u1 + dt*mu_t*(dtu1-(gamma_t*dtu));
				}
				B = [b, B[0]];
				Z = [z, Z[0]];
				dZ = [dz, dZ[0]];
				ddZ = [ddz, ddZ[0]];
			}
			u1
		};
		let ref dtu1 = self.dt(P, &u1);
		let ũ = map!(u1, &u, &dtu, dtu1; |u1, u, dtu, dtu1| u1-dt*(dtu+dtu1)/2.-u);
		let error = 0.8 * error(zip!(ũ, &u, u1).map(|(ũ, u, u1):(f64,&f64,&f64)| ũ / (absolute_tolerance + relative_tolerance*u.abs().max(u1.abs()))));
		if error > 1. {
			dt *= 0.8 / error.powf(1./3.);
			(v, jacobian_spectral_radius) = self.power_iteration(P, tmax, &u, &dtu, &v);
			//println!("{:e} {:e} {:e} {:e} {} {} {}", jacobian_spectral_radius, t, dt, error, nstep, stages, previous_error/error);
		} else {
			t += dt;
			if t >= tmax { break *u1; }
			u = *u1;
			dtu = *dtu1;
			nstep += 1;
			if nstep%25 == 0 { (v, jacobian_spectral_radius) = self.power_iteration(P, tmax, &u, &dtu, &v); }
			let factor = (0.8 * if previous_error > f64::EPSILON { dt/previous_dt*(previous_error/error).powf(1./3.) } else { 1./error.powf(1./3.) } ).clamp(0.1, 10.);
			if nstep%10000==0 { println!("{:e} {:e} {:e} {:e} {} {} {} {}", jacobian_spectral_radius, t, dt, error, nstep, stages, factor, previous_error/error); }
			previous_error = error;
			previous_dt = dt;
			dt = (factor*dt).min(tmax);
		}
	}
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
