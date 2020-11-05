#![allow(incomplete_features,non_snake_case,confusable_idents)]
#![feature(in_band_lifetimes,trait_alias,box_syntax,map_into_keys_values,associated_type_bounds,bindings_after_at,non_ascii_idents,min_const_generics,array_map,clamp)] //const_generics

pub fn scale(s: f64, v: impl IntoIterator<Item=f64,IntoIter:'t>) -> impl Iterator<Item=f64>+'t { v.into_iter().map(move |v| s*v) }
pub fn recip(x: impl IntoIterator<Item=f64,IntoIter:'t>) -> impl Iterator<Item=f64>+'t { x.into_iter().map(|x| f64::recip(x)) }
pub fn mul(a: impl IntoIterator<Item=f64,IntoIter:'t>, b: impl IntoIterator<Item=f64,IntoIter:'t>) -> impl Iterator<Item=f64>+'t { a.into_iter().zip(b.into_iter()).map(|(a,b)| a*b) }
pub fn dot(a: &[f64], b: &[f64]) -> f64 { mul(a.iter().copied(), b.iter().copied()).sum() }
pub fn dotia(a: impl IntoIterator<Item=f64>, b: &[f64]) -> f64 { mul(a, b.iter().copied()).sum() }
pub fn dotai(a: &[f64], b: impl IntoIterator<Item=f64>) -> f64 { mul(a.iter().copied(), b).sum() }

fn norm(iter: impl IntoIterator<Item=f64>) -> f64 { iter.into_iter().map(num::sq).sum::<f64>().sqrt() }
fn error(iter: impl ExactSizeIterator<Item=f64>) -> f64 { let len = iter.len(); (iter.sum::<f64>() / len as f64).sqrt() }

use itertools::izip;
//macro_rules! map { ($($args:ident),*; $expr:expr) => { itertools::izip!($($args,)*).map(|($($args,)*)| $expr) } }
//macro_rules! eval { ($($args:ident),*; $expr:expr) => { iter::array::FromIterator::from_iter(map!($($args),*; $expr)) } }
macro_rules! eval { ($($args:expr),*; $expr:expr) => { iter::array::FromIterator::from_iter(itertools::izip!($($args,)*).map($expr)) } }

mod ron {
use serde::Deserialize;
pub use std::collections::BTreeMap as Map;
#[derive(Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)] pub enum Element { H, O, Ar }
#[derive(Deserialize, Debug)] pub struct InitialState<'t> { pub temperature: f64, pub pressure: f64, #[serde(borrow)] pub mole_proportions: Map<&'t str, f64>, pub volume: f64 }
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
	#[serde(rename="beta")] pub temperature_exponent: f64,
	#[serde(rename="Ea")] pub activation_energy: f64
}

#[derive(Deserialize, Debug)] pub struct Troe { pub A: f64, pub T3: f64, pub T1: f64, pub T2: f64 }

#[derive(Deserialize, Debug)] pub enum Model<'t> {
	Elementary,
	ThreeBody { #[serde(borrow)] efficiencies: Map<&'t str, f64> },
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
	/*fn dT_Cp(&self, T: f64) -> f64 {
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		ideal_gas_constant * (a[1]+2.*a[2]*T+3.*a[3]*T*T+4.*a[4]*T*T*T)
	}
	fn dT_b(&self, T: f64) -> f64 { // dT(S/R - H/RT)
		use itertools::Itertools;
		let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
		(a[0]-1.)/T + a[1]/2. + a[2]/12.*T + a[3]/36.*T*T + a[4]/80.*T*T*T + a[5]/(T*T)
	}*/
}

pub fn arrhenius(&RateConstant{preexponential_factor, temperature_exponent, activation_energy}: &RateConstant, temperature: f64) -> f64 {
	preexponential_factor*temperature.powf(temperature_exponent)*f64::exp(-activation_energy/(ideal_gas_constant/4.184*temperature))
}

#[derive(Debug)] pub enum Model<const S: usize> {
	Elementary,
	ThreeBody { efficiencies: [f64; S] },
	Falloff { efficiencies: [f64; S], k0: RateConstant, troe: Troe },
}

impl<const S: usize> Model<S> {
fn efficiency(&self, T: f64, concentrations: &[f64; S], k_inf: f64) -> f64 {
	match self {
		Self::Elementary => 1.,
		Self::ThreeBody{efficiencies} => dot(efficiencies, concentrations),
		Self::Falloff{efficiencies, k0, troe: Troe{A, T3, T1, T2}} => {
			let Pr = dot(efficiencies, concentrations) * arrhenius(k0, T) / k_inf;
			let Fcent = (1.-A)*f64::exp(-T/T3)+A*f64::exp(-T/T1)+f64::exp(-T2/T);
			let log10Fcent = f64::log10(Fcent);
			let C = -0.4-0.67*log10Fcent;
			/*let N = 0.75-1.27*log10Fcent;
			let f1 = (f64::log10(Pr) + C)/(N-0.14*(f64::log10(Pr) + C));
			let F = num::exp10(log10Fcent/(1.+f1*f1));*/
			let ATroe = f64::log10(Pr) + C;
			let BTroe = 0.806 - 1.1762*log10Fcent - 0.14*f64::log10(Pr);
			let F = f64::powf(Fcent, f64::recip(1.+num::sq(ATroe/BTroe)));
			Pr * F / (1.+Pr) // Chemically activated bimolecular reaction
		}
	}
}
}

pub const S: usize = 9;
pub struct Reaction/*<const S: usize> where S>1*/ {
	pub equation: [(Box<[usize]>, Box<[u8]>); 2],
	rate_constant: RateConstant,
	pub model: Model<S>,
	specie_net_coefficients: [f64; S-1], // constant expression depends on a generic parameter
	//Σνf: f64,
	//Σνr: f64,
	PRν: f64,
}

pub struct System/*<const S: usize>*/ {
	pub molar_masses: [f64; S],
	pub thermodynamics: [NASA7; S],
	pub reactions: Box<[Reaction/*<S>*/]>,
}

macro_rules! Vec { () => ([f64; Self::N]) }
impl/*<const S: usize>*/ System/*<S> where S>1*/ {
	const N: usize = 2+S-1;
	//type Vec = [f64; Self::N]; // associated types are not yet supported in inherent impls (see #8995)

fn dt(&self, P: f64, y: &Vec!()) -> Vec!() {
	use iter::array::{from_iter, map, FromIterator};
	let Self{thermodynamics: species, reactions, molar_masses: W, ..} = self;

	let (T, V, n) = (y[0], y[1], &y[2..]);
	let T = T.clamp(100., 10000.0);
	let B = map(species, |s| s.b(T));
	let V = V.max(0.);
	let recipV = 1. / V;
	// C = n / V
	let concentrations = <[_;S-1]>::from_iter(n.iter().map(|n| recipV * n.max(0.)));
	let C = P / (ideal_gas_constant * T);
	let concentrations = from_iter(concentrations.iter().copied().chain(std::iter::once(C - concentrations.iter().sum::<f64>())));

	let a = S-1;
	let mut dtω = vec!(0.; S-1); // Skips most abundant specie (last index) (will be deduced from conservation)
	for Reaction{equation, rate_constant, model, specie_net_coefficients: ν, PRν, ..} in reactions.iter() {
		let equilibrium_constant = PRν * f64::exp(dot(ν, &B));
		let kf = arrhenius(rate_constant, T);
		let kr = kf / equilibrium_constant;
		// ΠC^ν
		let [ΠCνf, ΠCνr] = from_iter(equation.iter().map(|(species, ν)| species.iter().zip(ν.iter()).map(|(&specie, &ν)| concentrations[specie].powi(ν as i32)).product::<f64>() ));
		let Rf = kf * ΠCνf;
		let Rr = kr * ΠCνr;
		let c = model.efficiency(T, &concentrations, kf);
		let R = Rf - Rr;
		let cR = c * R;
		for (specie, ν) in ν.iter().enumerate() {
			//ω̇̇̇̇̇̇̇̇̇̇ = Σ ν c R
			dtω[specie] += ν * cR;
		}
	}
	let Cp = map(species, |s| s.specific_heat_capacity(T));
	// 1 / (Σ C.Cp)
	let rcp_ΣCCp = 1./dot(&concentrations, &Cp);
	let H = map(species, |s| s.specific_enthalpy(T));
	// dtT = - 1 / (Σ C.Cp) . Σ H.ω̇
	let dtT = - rcp_ΣCCp * dot(&H, &dtω);
	// dt(V) = V.(dt(T) / T + TR/P.Σ{(1-Wk/Wa).ω})
	let Wa = W[a];
	let dtE = dotia(W.iter().map(|W| 1.-W/Wa), &dtω);
	let dtV = V * (dtT / T + T * ideal_gas_constant / P * dtE);
	// dt(n) = Vω
	let dtn = dtω.iter().map(|dtω| V*dtω);
	from_iter([dtT, dtV].iter().copied().chain(dtn))
}

fn spectral_radius(&self, P: f64, tmax: f64, y: &Vec!(), dty: &Vec!(), v: &mut Vec!()) -> f64 {
	let [norm_y, norm_v] = [y,v].map(|x| norm(x.iter().copied()));
	// ?
	let dynrm = match [norm_y != 0., norm_v != 0.] {
		[true, true] => {
			let dynrm = norm_y * f64::EPSILON.sqrt();
			for (v, y) in izip!(&mut *v, y) { *v = y + *v * dynrm / norm_v; }
			dynrm
		}
		[true, false] => {
			for (v, y) in izip!(&mut *v, y) { *v = y * (1. + f64::EPSILON.sqrt()); }
			norm_y * f64::EPSILON.sqrt()
		}
		[false, true] => {
			let dynrm = f64::EPSILON;
			for v in &mut *v  { *v *= dynrm / norm_v; }
			dynrm
		}
		[false, false] => {
			for v in &mut *v { *v = f64::EPSILON }
			f64::EPSILON
		}
	};
	// nonlinear power method
	let mut sigma = 0.;
	for i in 1..=50 {
		let dtv = self.dt(P, v);
		let norm = norm(izip!(&dtv, dty).map(|(dtv, dty)| dtv - dty));
		let previous_sigma = sigma;
		sigma = norm / dynrm;
		if i >= 2 && (sigma - previous_sigma).abs() <= 0.01*sigma.max(1./tmax) {
			for (v, y) in izip!(&mut *v, y) { *v = *v - y; }
			break;
		}
		if norm != 0. {
			for (v, y, dtv, dty) in izip!(&mut *v, y, &dtv, dty) { *v = y + (dtv - dty) * (dynrm / norm); }
		} else {
				let i = i % Self::N;
				v[i] = y[i] - (v[i] - y[i]);
				panic!();
		}
	}
	sigma * 1.2
}

pub fn integrate(&self, rtol: f64, atol: f64, tmax: f64, P: f64, y: &mut Vec!()) {
	let max_steps = ((rtol / (10. * f64::EPSILON)).sqrt().round() as usize).max(2);
	let mut nstep = 0;
	let mut t = 0.;
	let mut dty = self.dt(P, y);
	let mut v = dty;
	let mut jacobian_spectral_radius = self.spectral_radius(P, tmax, y, &dty, &mut v);
	let mut dt = {
		let dt = (1./jacobian_spectral_radius).min(tmax);
		let dty1 = self.dt(P, &eval!(&*y, &dty; |(y, dty)| *y + (dt * dty)));
		(dt/(dt*error(izip!(&dty, &dty1, &*y).map(|(dty, dty1, y)| (dty1 - dty) / (atol + rtol * y.abs())))) / 10.).min(tmax)
	};
	let (mut previous_error, mut previous_dt) = (0., 0.);
	loop {
		if 1.1*dt >= tmax - t { dt = tmax- t; } // fit last step
		let steps = 1 + (1. + 1.54 * dt * jacobian_spectral_radius).sqrt().floor() as usize;
		let steps = if steps > max_steps {
			dt = (max_steps*max_steps - 1) as f64 / (1.54 * jacobian_spectral_radius);
			max_steps
		} else { steps };
		let w0 = 1. + 2. / (13.0 * (steps * steps) as f64);
		let sqw01 = w0*w0 - 1.;
		let arg = steps as f64 * (w0 + sqw01.sqrt()).ln();
		let w1 = arg.sinh() * sqw01 / (arg.cosh() * steps as f64 * sqw01.sqrt() - w0 * arg.sinh());
		let mut B = [1. / (4.*w0*w0); 2];
		let mu_t = w1 * B[0];
		let [mut y0, mut y1] = [*y, eval!(&*y, &dty; |(y, dty)| y + mu_t * dt * dty)];
		let mut Z = [w0, 1.];
		let mut dZ = [1., 0.];
		let mut ddZ = [0., 0.];
		let mut steps = steps - 2;
		loop {
			let z = 2. * w0 * Z[0] - Z[1];
			let dz = 2. * w0 * dZ[0] - dZ[1] + 2. * Z[0];
			let ddz = 2. * w0 * ddZ[0] - ddZ[1] + 4. * dZ[0];
			let b = ddz / (dz * dz);
			let gamma_t = 1. - (Z[0] * B[0]);
			let nu = - b / B[1];
			let mu = 2. * b * w0 / B[0];
			let mu_t = mu * w1 / w0;
			let dty1 = self.dt(P, &y1);
			for (y0, y1, dty1, y, dty) in izip!(&mut y0, &mut y1, &dty1, &*y, &dty) {
				let y0_ = *y0;
				*y0 = *y1;
				*y1 = (1.-mu-nu)*y + nu*y0_ + mu**y1 + dt*mu_t*(dty1-(gamma_t*dty));
			}
			if steps == 0 { break; }
			steps -= 1;
			B = [b, B[0]];
			Z = [z, Z[0]];
			dZ = [dz, dZ[0]];
			ddZ = [ddz, ddZ[0]];
		}
		let dty1 = self.dt(P, &y1);
		let error = error(izip!(&*y,&y1,&dty,&dty1).map(|(y,y1,dty,dty1)| (0.8*(y1-y)+0.4*dt*(dty+dty1)/(atol + rtol*y.abs().max(y1.abs()))).powi(2)));
		if error > 1. { // error too large, step is rejected
			dt *= 0.8 / error.powf(1./3.);
			assert!(dt >= f64::EPSILON);
			jacobian_spectral_radius = self.spectral_radius(P, tmax, y, &dty, &mut v);
		} else { // step accepted
			t += dt;
			*y = y1;
			if t >= tmax { break; }
			dty = dty1;
			nstep += 1;
			if (nstep % 25) == 0 { jacobian_spectral_radius = self.spectral_radius(P, tmax, y, &dty, &mut v); }
			let factor = (0.8 * if previous_error > f64::EPSILON { dt/previous_dt*previous_error.powf(1./3.) } else { 1. } / error.powf(1./3.)).clamp(0.1, 10.);
			previous_error = error;
			previous_dt = dt;
			dt *= factor;
		}
	}
}
}

/*impl<const S: usize> State<S> {
	pub fn step(&mut self, System{thermodynamics: species, reactions, amount, volume: V, molar_masses: W}: &System<S>) {
		use iter::array::{from_iter, map, generate};
		let scale = |s, v| from_iter(scale(s, v.iter().copied()));
		macro_rules! map { ($($args:ident),*| $expr:expr) => { izip!($($args),*).map(|($($args),*)| $expr) } }

		let rcpn = 1. / amount;
		let T = self.temperature;
		let B = map(species, |s| s.b(T));
		let dTB = map(species, |s| s.dT_b(T));
		let rcpV = 1. / V;
		// C = n / V
		let C = scale(rcpV, self.amounts);
		let a = S-1;
		let Ca = C[a];
		let mut ω = vec!(0.; S-1); // Skips most abundant specie (last index) (will be deduced from conservation)
		let mut dTω = vec!(0.; S-1);
		let mut dVω = vec!(0.; S-1);
		let mut dnω = vec!(vec!(0.; S-1); S-1); //[k][j]
		for Reaction{equation, rate_constant: rate_constant@RateConstant{temperature_exponent: β, activation_energy: Ea}, model, specie_net_coefficients: ν, Σνf, Σνr, PRν} in reactions.iter() {
			let equilibrium_constant = PRν * f64::exp(dot(ν, B));
			let kf = arrhenius(rate_constant, T);
			let kr = kf / equilibrium_constant;
			// ΠC^ν
			let [ΠCνf, ΠCνr] = from_iter(equation.iter().map(|(species, ν)| species.iter().zip(ν.iter()).map(|(&specie, &ν)| C[specie].powi(ν as i32)).product::<f64>() ));
			let Rf = kf * ΠCνf;
			let Rr = kr * ΠCνr;
			let νfRfνrRr = vec!(0.; S);
			let [forward, reverse] = equation;
			for (specie, ν) in forward { νfRfνrRr[specie] += ν * Rf; }
			for (specie, ν) in reverse { νfRfνrRr[specie] -= ν * Rr; }
			let c = model.efficiency(T, &C, kf);
			let R = Rf - Rr;
			// dT(R) = (β+Ea/(T.Rideal))/T.R + Rr. Σ ν.dT(B) - νfRfνrRr[a]/T
			let dTR = (β+Ea/(T*ideal_gas_constant))/T*R + Rr*dot(ν,dTB) - νfRfνrRr[a] / T;
			// dV(R) = 1/V . ( (kf.Sf-kr.Sr) - (Σνf.Rf - Σνr.Rr) )
			let dVR = rcpV * ( νfRfνrRr[a] - (Σνf*Rf - Σνr*Rr));
			// dn(R) = 1/n . ( kf.(Sf-Sfa) - kr.(Sr-Sra) )
			let dnR = map(νfRfνrRr, |νfRfνrRrj| rcpn * (νfRfνrRrj - νfRfνrRr[a]));
			let (dTc, dVc, dnc) = match model {
				Model::Elementary => (0., 0., vec!(0.; S-1)),
				Model::ThirdBody{efficiencies}|Model::Falloff{efficiencies} => {
					let has = map(efficiencies, |e| if e != 0. { 1. } else { 0. });
					(
						// dT(c) = has(a) . (-C/T)
						has[a] * -C/T,
						// dV(c) = 1/V . ( has(a). C - Σ C )
						rcpV * (has[a] * C  - dot(has, C)),
						// dn(c) = 1/V . ( has(n) - has(a) )
						has[..a].map(|has_n| rcpV * (has_n - has[a]))
					)
				}
			};
			let cR = c * R;
			let RdTccdTR = R * dTc + c * dTR;
			let RdVccdVR = R * dVc + c * dVR;
			let RdnccdnR = from_iter(map!(dnc,dnR| R*dnc + c*dnR));
			for (specie, ν) in ν.iter().enumerate() {
				//ω̇̇̇̇̇̇̇̇̇̇ = Σ ν c R
				ω[specie] += ν * cR;
				// dT(ω̇̇̇̇̇̇̇̇̇̇) = Σ ν.(R.dT(c)+c.dT(R))
				dTω[specie] += ν * RdTccdTR;
				// dV(ω) = Σ ν.(R.dV(c)+c.dV(R))
				dVω[specie] += ν * RdVccdVR;
				// dn(ω) = Σ ν.(R.dn(c)+c.dn(R))
				for dnω in dnω[specie] { dnω += ν * RdnccdnR; }
			}
		}
		use nalgebra::{Matrix, MatrixMN};
		let mut J = unsafe{MatrixMN::<f64, {2+S-1}, {2+S-1}>::new_uninitialized()}; // fixme
		let Cp = map(species, |s| s.specific_heat_capacity(T));
		// 1 / (Σ C.Cp)
		let rcp_ΣCCp = 1./dot(C, Cp);
		let H = species.map(|s| s.specific_enthalpy(T));
		let Wa = W[a];
		let HaWa = H[a]/Wa;
		// Ha/Wa*W - H
		let HaWaWH = from_iter(map!(W,H| HaWa*W - H));
		// dtT = - 1 / (Σ C.Cp) . Σ H.ω̇
		let dtT = - rcp_ΣCCp * dot(H, ω);
		let CpaWa = Cp[a]/Wa;
		// dT(dtT) = 1 / (Σ C.Cp) . [ dtT . Σ C.(Cpa/T - dT(Cp)) + Σ_ ( (Ha/Wa*W - H).dT(ω) + (Cpa/Wa.W - Cp).ω ) ]
		let dTdtT = rcp_ΣCCp * ( dtT * dot(C, species.iter().map(|s| Cp[a]/T - s.dt_Cp())) + dot(HaWaWH, dTω) + dotia(map!(W,Cp, CpaWa*W - Cp), ω));
		J[(0,0)] = dTdtT;
		// dV(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).dV(ω) + dtT/V . Σ_ C.(Cp-Cpa) ]
		let dVdtT = rcp_ΣCCp * ( dot(HaWaWH, dVω) + rcpV * dtT * dotai(C, Cp.iter().map(|Cp| Cp - Cp[a])));
		J[(0,1)] = dVdtT;
		// dn(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).dn(ω) + dtT/V . (Cpa-Cp) ]
		let dndtT = from_iter(map!(dnω, Cp| rcp_ΣCCp * ( dot(HaWaWH, dnω) + rcpV * dtT * (Cp[a]-Cp) )));
		J.row_part_mut(0,2,S+1).copy_from_slice(dndtT);

		// dT(dtV) = V/C . Σ_ (1-W/Wa).(dT(ω)+ω/T) + V/T.(dT(dtT) - dtT/T)
		let rcpC = V*rcpn;
		let WWa = map(W, |W| 1.-W/Wa);
		let VT = V/T;
		let dTdtV = V*rcpC * map!(WWa,dTω,ω| WWa*(dTω+ω/T)).sum() + VT * (dTdtT - dtT/T);
		J[(1,0)] = dTdtV;
		// dV(dtn) = VdV(ω)+ω
		let dVdtn = from_iter(map!(dVω,ω| V*dVω+ω));
		// dV(dtV) = 1/C . Σ_ (1-W/Wa).dV(dtn) + 1/T.(V.dV(dtT)+dtT)
		let dVdtV = rcpC * dot(WWa, dVdtn) + 1./T*(V*dVdtT+dtT);
		J[(1,1)] = dVdtV;
		// dn(dtn) = Vdn(ω)
		let dndtn = generate(S-1, |j| generate(S-1, |k| V*dnω[k][j])); // Transpose [k][j] -> [j][k]
		// dn(dtV) = 1/C . Σ_ (1-W/Wa).dn(dtn) + V/T.dn(dtT))
		let dndtV = map!(dndtn, dndtT| rcpC * dot(WWa, dndtn)+ VT*dndtT);
		J.row_part_mut(1,2,S+1) = Matrix::from_iterator(dndtV);

		// dT(dtn) = VdT(ω)
		let dTdtn = scale(V, dTω);
		J.column_part_mut(0,2,S+1) = Matrix::from_iterator(dTdtn);
		// dV(dtn)
		J.column_part_mut(1,2,S+1) = Matrix::from_iterator(dVdtn);
		// dn(dtn)
		for j in 2..S+1 { J.column_part_mut(j,2,S+1).copy_from_slice(dndtn[j]); }

		// dt(V) = V.(TR/P.Σ{(1-Wk/Wa).ω}+1/T.dt(T))
		// dt(T) = - Σ_ ((H-(W/Wa).Ha).w) / ( Ct.Cpa + Σ_ (Cp_Cpa)C )
		// dt(n) = Vω
		for (state, rate) in self.amounts[..specie_count-1].iter_mut().zip(production_rates.iter()) { *state = 0f64.max(*state + system.time_step * system.volume * rate); }
		let total_amount = system.pressure * system.volume / (ideal_gas_constant * self.temperature);
		self.amounts[specie_count-1] = 0f64.max(total_amount - self.amounts[..specie_count-1].iter().sum::<f64>()); // Enforces mass conservation constraint by rescaling most abundant specie (last index)
		for &amount in self.amounts.iter() { assert!(amount>= 0. && amount< total_amount,"{} {:?} {}", amount, &self.amounts, total_amount); }
		self.temperature += system.time_step * dtT;
	}
}*/

#[derive(Clone)] pub struct State<const S: usize> {
	pub temperature: f64,
	//pub volume: f64,
	pub amounts: [f64; S]
}
impl<const S: usize> State<S> {
	pub const S: usize = S;
}

pub struct Simulation<'t, const S: usize> {
	pub species: Box<[&'t str]>,
	pub system: System/*<S>*/,
	pub amount: f64,
	pub time_step: f64,
	pub pressure: f64,
	pub volume: f64,
	pub state: State<S>
}

impl<const S: usize> Simulation<'t, S> {
#[fehler::throws(anyhow::Error)] pub fn new(system: &'b [u8]) -> Self where 'b: 't {
	let ron::System{species, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;

	let standard_atomic_weights : Map<Element, f64> = ::ron::de::from_str("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}")?;
	let standard_atomic_weights : Map<_, f64> = standard_atomic_weights.into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect();

	use iter::{from_iter, array::{self, map, generate}};
	let specie_names = iter::from_iter(species.keys().copied());
	let molar_masses = array::from_iter(species.values().map(|Specie{composition,..}| composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum()));
	let thermodynamics = array::from_iter(species.into_values().map(|Specie{thermodynamic,..}| thermodynamic));
	let species = specie_names;

	let reactions = iter::from_iter(Vec::from(reactions).into_iter().map(|self::ron::Reaction{equation, rate_constant, model}| {
		let atmospheric_pressure = 101325.;
		//let equation = iter::array::Iterator::collect::<[(Box<[_]>, Box<[_]>);2]>(equation.iter().map(|e| (e.keys().map(|&key| species.iter().position(|&k| k==key).expect(key)).collect(), e.values().copied().collect())));
		let equation = map(&equation, |e| (iter::from_iter(e.keys().map(|&key| species.iter().position(|&k| k==key).expect(key))), iter::from_iter(e.values().copied())));
		let specie_net_coefficients = generate(|specie| {let [reactant, product] = map(&equation, |(species, ν)| species.iter().position(|&s| s==specie).map(|i| ν[i] as i8).unwrap_or(0)); (product-reactant) as f64});
		//let [(_, νf), (_, νr)] = &equation;
		//let [Σνf, Σνr] = [νf.iter().sum::<u8>() as f64, νr.iter().sum::<u8>() as f64];
		let PRν = f64::powf(atmospheric_pressure / ideal_gas_constant, specie_net_coefficients.iter().sum());
		Reaction{
			equation,
			rate_constant,
			model: {use self::ron::Model::*; match model {
				Elementary => Model::Elementary,
				ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: array::from_iter(species.iter().map(|specie| *efficiencies.get(specie).unwrap_or(&1.)))},
				Falloff{efficiencies, k0, troe} => Model::Falloff{efficiencies: array::from_iter(species.iter().map(|specie| *efficiencies.get(specie).unwrap_or(&1.))), k0, troe},
			}},
			specie_net_coefficients,
			//Σνf, Σνr,
			PRν
		}
	}));

	let Phase::IdealGas{state, ..} = Vec::from(phases).into_iter().next().unwrap();
	let InitialState{temperature, pressure, mole_proportions, volume} = state;
	let mole_proportions = from_iter(species.iter().map(|specie| *mole_proportions.get(specie).unwrap_or(&0.)));
	let amount = pressure * volume / (ideal_gas_constant * temperature);
	let amounts = array::from_iter(scale(amount/mole_proportions.iter().sum::<f64>(), mole_proportions.iter().copied()));

	Self{
		species,
		system: System{molar_masses, thermodynamics, reactions},
		amount, time_step, pressure, volume,
		state: State{temperature, /*volume,*/ amounts}
	}
}
}

impl<const S: usize> From<State<S>> for Box<[Box<[f64]>]> { fn from(s: State<S>) -> Self { box [box [s.temperature] as Box<[_]>, box s.amounts] as Box<[_]> } }
