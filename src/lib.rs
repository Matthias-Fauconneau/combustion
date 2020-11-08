#![allow(incomplete_features,non_snake_case,confusable_idents)]
#![feature(in_band_lifetimes,trait_alias,box_syntax,map_into_keys_values,associated_type_bounds,bindings_after_at,non_ascii_idents,min_const_generics,clamp,array_value_iter,type_ascription,array_map,once_cell)]
#![no_implicit_prelude]
#![allow(non_upper_case_globals)]
//mod array;
mod ron;
pub const ideal_gas_constant : f64 = 8.31446261815324;
extern crate fehler; extern crate anyhow;
extern crate iter; use iter::{map, eval, zip, Zip, IntoChain, Dot, IntoEnumerate, Sub, generate, collect, r#box};
extern crate std; use std::{convert::{From, TryInto}, iter::{IntoIterator, Iterator}, panic, boxed::Box, cmp::Ord, matches};
extern crate num; use num::norm;
pub fn ssq(iter: impl iter::IntoExactSizeIterator+IntoIterator<Item=f64>) -> f64 {
	let iter = iter.into_iter();
	let len = std::iter::ExactSizeIterator::len(&iter);
	f64::sqrt(std::iter::Iterator::sum::<f64>(std::iter::Iterator::map(iter, num::sq)) / len as f64)
}
pub trait Suffix<T> { fn suffix<const S: usize>(&self) -> &[T; S]; }
impl<T, const N: usize> Suffix<T> for [T; N] { fn suffix<const S: usize>(&self) -> &[T; S] { std::convert::TryFrom::try_from(&self[N-S..]).unwrap() } }
//use array::{Suffix}; //Array;
use self::ron::*;

pub struct NASA7([[f64; 7]; 2]);
impl NASA7 {
	const T_split : f64 = 1000.;
	pub fn a(&self, T: f64) -> &[f64; 7] { if T < Self::T_split { &self.0[0] } else { &self.0[1] } }
	pub fn specific_heat_capacity(&self, T: f64) -> f64 { let a = self.a(T); ideal_gas_constant * (a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T) }
	pub fn specific_enthalpy(&self, T: f64) -> f64 { let a = self.a(T); ideal_gas_constant * (a[5]+a[0]*T+a[1]/2.*T*T+a[2]/3.*T*T*T+a[3]/4.*T*T*T*T+a[4]/5.*T*T*T*T*T) }
	pub fn specific_entropy(&self, T: f64) -> f64 { let a = self.a(T); ideal_gas_constant * (a[6]+a[0]*f64::ln(T)+a[1]*T+a[2]/2.*T*T+a[3]/3.*T*T*T+a[4]/4.*T*T*T*T) }
	pub fn b(&self, T: f64) -> f64 { let a = self.a(T); a[6] - a[0] + (a[0]-1.)*f64::ln(T) + a[1]/2.*T + a[2]/6.*T*T + a[3]/12.*T*T*T + a[4]/20.*T*T*T*T - a[5]/T } // S/R - H/RT
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
		Self::ThreeBody{efficiencies} => efficiencies.dot(concentrations),
		Self::Falloff{efficiencies, k0, troe: Troe{A, T3, T1, T2}} => {
			let Pr = efficiencies.dot(concentrations) * arrhenius(k0, T) / k_inf;
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

pub struct Reaction<const S: usize> {
	pub equation: [Zip<Box<[usize]>, Box<[u8]>>; 2],
	rate_constant: RateConstant,
	pub model: Model<S>,
	specie_net_coefficients: [f64; S],
	//Σνf: f64,
	//Σνr: f64,
	PRν: f64,
}

pub struct System<const S: usize> {
	pub molar_masses: [f64; S],
	pub thermodynamics: [NASA7; S],
	pub reactions: Box<[Reaction<S>]>,
}

impl<const S: usize> System<S> {
fn dt(&self, P: f64, y: &[f64; S]) -> [f64; S] {
	let Self{thermodynamics: species, reactions, molar_masses: W, ..} = self;

	let (T, V, n) = (y[0], y[1], y.suffix());
	let T = T.clamp(100., 10000.0);
	let B = eval!(species; |s| s.b(T));
	let V = V.max(0.);
	let recipV = 1. / V;
	let concentrations : [_; /*S-1*/8] = eval!(n; |n| recipV * n.max(0.));
	let C = P / (ideal_gas_constant * T);
	let concentrations : [_; S] = collect(concentrations.chain([C - concentrations.iter().sum():f64]));

	let a = S-1;
	let ref mut dtω = [0.; S][..S-1]; // Skips most abundant specie (last index) (will be deduced from conservation)
	for Reaction{equation, rate_constant, model, specie_net_coefficients: ν, PRν, ..} in reactions.iter() {
		let equilibrium_constant = PRν * f64::exp(ν.dot(&B));
		let kf = arrhenius(rate_constant, T);
		let kr = kf / equilibrium_constant;
		let [ΠCνf, ΠCνr] : [f64;2] = eval!(equation; |side| side.a.into_iter().zip(side.b.into_iter()).map(|(specie, ν)| concentrations[*specie].powi(*ν as i32)).product());
		let Rf = kf * ΠCνf;
		let Rr = kr * ΠCνr;
		let c = model.efficiency(T, &concentrations, kf);
		let R = Rf - Rr;
		let cR = c * R;
		for (specie, ν) in ν.enumerate() { dtω[specie] += ν * cR; }
	}
	let ref dtω = *dtω;
	let Cp = map!(species; |s| s.specific_heat_capacity(T));
	let rcp_ΣCCp = 1./concentrations.dot(Cp);
	let H = map!(species; |s| s.specific_enthalpy(T));
	let dtT = - rcp_ΣCCp * dtω.dot(H);
	let Wa = W[a];
	let dtE = W.map(|W| 1.-W/Wa).dot(dtω);
	let dtV = V * (dtT / T + T * ideal_gas_constant / P * dtE);
	let dtn = map!(dtω; |dtω| V*dtω);
	collect([dtT, dtV].chain(dtn))
}

// Estimate principal eigenvector/value of dyF|y
fn power_iteration(&self, P: f64, tmax: f64, y: &[f64; S], dty: &[f64; S], v: &[f64; S]) -> ([f64; S], f64) {
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
		if i >= 2 && f64::abs(ρ - previous_ρ) <= 0.01*ρ.max(1./tmax) { break; } // Early exit
		*yεv = eval!(y, dtεv; |y, dtεv| y + (ε / norm_dtεv) * dtεv);
		//yεv = y.zip(&dtεv).map(|(y, dtεv)| y + (ε / norm_dtεv) * dtεv);
	}
	(yεv.sub(y), ρ * 1.2)
}

pub fn step(&self, rtol: f64, atol: f64, tmax: f64, P: f64, mut y: [f64; S]) -> [f64; S] {
	let max_steps = ((rtol / (10. * f64::EPSILON)).sqrt().round() as usize).max(2);
	let mut nstep = 0;
	let mut t = 0.;
	let ref mut dty = self.dt(P, &y);
	let (mut v, mut jacobian_spectral_radius) = self.power_iteration(P, tmax, &y, dty, dty);
	let mut dt = {
		let dt = (1./jacobian_spectral_radius).min(tmax);
		let ref dty1 = self.dt(P, &eval!(&y, &*dty; |y, dty| y + dt * dty));
		(dt/(dt*ssq(map!(&*dty, dty1, &y; |dty, dty1, y| (dty1 - dty) / (atol + rtol * y.abs())))) / 10.).min(tmax)
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
		let [ref mut y0, mut y1] = [y, eval!(&y, &*dty; |y, dty| y + mu_t * dt * dty)];
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
			let ref dty1 = self.dt(P, &y1);
			for (y0, y1, dty1, y, dty) in zip!(y0, &mut y1, dty1, &y, &*dty) {
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
		let ref dty1 = self.dt(P, &y1);
		let error = ssq(map!(&y,&y1,&*dty,dty1; |y,y1,dty,dty1| (0.8*(y1-y)+0.4*dt*(dty+dty1)/(atol + rtol*y.abs().max(y1.abs()))).powi(2)));
		if error > 1. { // error too large, step is rejected
			dt *= 0.8 / error.powf(1./3.);
			assert!(dt >= f64::EPSILON);
			{let t = self.power_iteration(P, tmax, &y, &dty, &v); v = t.0; jacobian_spectral_radius = t.1;}
		} else { // step accepted
			t += dt;
			if t >= tmax { break y1; }
			y = y1;
			*dty = *dty1;
			nstep += 1;
			if (nstep % 25) == 0 {let t = self.power_iteration(P, tmax, &y, &dty, &v); v = t.0; jacobian_spectral_radius = t.1;}
			let factor = (0.8 * if previous_error > f64::EPSILON { dt/previous_dt*previous_error.powf(1./3.) } else { 1. } / error.powf(1./3.)).clamp(0.1, 10.);
			previous_error = error;
			previous_dt = dt;
			dt *= factor;
		}
	}
}
}

#[derive(Clone)] pub struct State<const S: usize> {
	pub temperature: f64,
	//pub volume: f64,
	pub amounts: [f64; S]
}
impl<const S: usize> State<S> {
	pub const S: usize = S;
}

pub struct Simulation<'t, const S: usize> {
	pub species: [&'t str; S],
	pub system: System<S>,
	pub amount: f64,
	pub time_step: f64,
	pub pressure: f64,
	pub volume: f64,
	pub state: State<S>
}

use std::lazy::SyncLazy;
static standard_atomic_weights : SyncLazy<Map<Element, f64>> = SyncLazy::new(|| {
	Iterator::collect((::ron::de::from_str("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}").unwrap():Map<Element, f64>).into_iter().map(|(e,g)| (e, g*1e-3/*kg/g*/)))
});

impl<const S: usize> Simulation<'t, S> {
#[fehler::throws(anyhow::Error)] pub fn new(system: &'b [u8]) -> Self where 'b: 't {
	let ron::System{species, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;

	let specie_names : [_; S] = collect(species.keys().copied());
	let molar_masses = collect(species.values().map(|Specie{composition,..}| composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum()));
	let thermodynamics = collect(species.into_values().map(|Specie{thermodynamic:ron::NASA7{temperature_ranges,coefficients},..}| {
		assert!(matches!(temperature_ranges[..], [_,Tsplit,_] if Tsplit == NASA7::T_split));
		NASA7(coefficients[..].try_into().unwrap())
	}));
	let species = specie_names;

	let reactions = r#box::collect(reactions.into_vec().into_iter().map(|self::ron::Reaction{ref equation, rate_constant, model}| {
		let atmospheric_pressure = 101325.;
		let equation = eval!(equation; |e| iter::IntoZip::zip(r#box::collect(e.keys().map(|&key| species.iter().position(|&k| k==key).expect(key))), r#box::collect(e.values().copied())));
		let specie_net_coefficients = generate(|specie| {
			let [reactant, product]:[_;2] = collect(equation.iter().map(|side| side.a.into_iter().zip(side.b.into_iter()).find(|(&s,_)| s==specie).map(|(_,&ν)| ν as i8).unwrap_or(0)));
			(product-reactant) as f64
		});
		let PRν = f64::powf(atmospheric_pressure / ideal_gas_constant, specie_net_coefficients.iter().sum());
		Reaction{
			equation,
			rate_constant,
			model: {use {iter::VectorCollect, self::ron::Model::*}; match model {
				Elementary => Model::Elementary,
				ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: iter::IntoMap::map(&species, |specie| *efficiencies.get(specie).unwrap_or(&1.)).collect()},
				Falloff{efficiencies, k0, troe} => Model::Falloff{efficiencies: iter::IntoMap::map(&species, |specie| *efficiencies.get(specie).unwrap_or(&1.)).collect(), k0, troe},
			}},
			specie_net_coefficients,
			//Σνf, Σνr,
			PRν
		}
	}));

	let Phase::IdealGas{state, ..} = phases.into_vec().into_iter().next().unwrap();
	let InitialState{temperature, pressure, mole_proportions, volume} = state;
	let mole_proportions = eval!(&species; |specie| *mole_proportions.get(specie).unwrap_or(&0.));
	let amount = pressure * volume / (ideal_gas_constant * temperature);
	let amounts = eval!(&mole_proportions; |mole_proportion| amount/mole_proportions.iter().sum::<f64>() * mole_proportion);

	Self{
		species,
		system: System{molar_masses, thermodynamics, reactions},
		amount, time_step, pressure, volume,
		state: State{temperature, /*volume,*/ amounts}
	}
}
}

impl<const S: usize> From<State<S>> for Box<[Box<[f64]>]> { fn from(s: State<S>) -> Self { box [box [s.temperature] as Box<[_]>, box s.amounts] as Box<[_]> } }
