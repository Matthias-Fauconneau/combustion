#![feature(min_const_generics,non_ascii_idents,in_band_lifetimes,once_cell,clamp,array_map,map_into_keys_values)]
#![allow(non_snake_case,confusable_idents,non_upper_case_globals)]
#[macro_export] macro_rules! assert { ($cond:expr $(, $val:expr)* ) => { std::assert!($cond, "{}. {:?}", stringify!($cond), ( $( format!("{} = {:?}", stringify!($val), $val), )* ) ); } }
pub mod ron;
use {std::convert::TryInto,
				iter::{IntoIterator, array_from_iter as from_iter, box_collect, into::{Collect, Enumerate, IntoChain, Zip, IntoZip, IntoMap, map, Find, Product}, zip, eval,
				vec::{eval, generate, Sub, Dot, Suffix}},
				num::{ssq, norm}};
trait IntoFormat : IntoIterator+Sized {
	fn format(self, sep: &str) -> itertools::Format<'_, Self::IntoIter> { itertools::Itertools::format(self.into_iter(), sep) }
	fn format_with<F: FnMut(Self::Item, &mut dyn FnMut(&dyn std::fmt::Display) -> std::fmt::Result) -> std::fmt::Result>(self, sep: &str, format: F) -> itertools::FormatWith<'_, Self::IntoIter, F> {
		itertools::Itertools::format_with(self.into_iter(), sep, format)
	}
}
impl<I:IntoIterator> IntoFormat for I {}

pub fn error<I:iter::IntoExactSizeIterator+iter::IntoIterator<Item=f64>>(iter: I) -> f64 {
	let iter = iter::IntoIterator::into_iter(iter);
	let len = iter.len();
	(ssq(iter) / len as f64).sqrt()
}

pub const ideal_gas_constant : f64 = 8.31446261815324;

pub use self::ron::*;

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

pub struct Reaction<const S: usize, const S1: usize> {
	pub equation: [Zip<Box<[usize]>, Box<[u8]>>; 2],
	rate_constant: RateConstant,
	pub model: Model<S>,
	specie_net_coefficients: [f64; S1],
	//Σνf: f64,
	//Σνr: f64,
	PRν: f64,
}

pub struct System<const S: usize, const S1: usize, const N: usize> {
	pub molar_masses: [f64; S],
	pub thermodynamics: [NASA7; S],
	pub reactions: Box<[Reaction<S,S1>]>,
}

//extern "C" { fn cantera(rtol: f64, atol: f64, T: f64, P: f64, X: *const std::os::raw::c_char, len: &mut usize, species: &mut *const *const std::os::raw::c_char, dtw: &mut *const f64); }

impl<const S: usize, const S1: usize, const N: usize> System<S,S1,N> {
fn dt(&self, P: f64, y: &[f64; N]) -> [f64; N] {
	//let specie_names = ["H","H2","O","OH","H2O","O2","HO2","H2O2","AR"]; let specie = |name:&str| specie_names.iter().position(|key|key==&name); use itertools::Itertools; // DEBUG
	let Self{thermodynamics: species, reactions, molar_masses: W, ..} = self;

	let (T, V, n) = (y[0], y[1], y.suffix());
	let recipV = 1. / V;
	let concentrations : [_; /*S-1*/S1] = eval(n, |n| recipV * n); // Skips most abundant specie (last index) (will be deduced from conservation)
	let C = P / (ideal_gas_constant * T);
	let concentrations : [_; S] = from_iter(concentrations.chain([C - concentrations.iter().sum::<f64>()]));

	let B = eval(species, |s| s.b(T));
	let ref mut dtω = [0.; /*S-1*/S1]; //][..S-1];
	for Reaction{equation, rate_constant, model, specie_net_coefficients: ν, PRν, ..} in reactions.iter() {
		let equilibrium_constant = PRν * f64::exp(ν.dot(&B));
		let kf = arrhenius(rate_constant, T);
		let kr = kf / equilibrium_constant;
		//println!("{}", equation.format_with(" = ", |side, f| f(&side.format_with(" + ", |(&specie, ν), f| f(&format_args!("{} {} {}", ν, specie_names[specie], concentrations[specie]))))));
		let [ΠCνf, ΠCνr] : [f64;2] = eval(equation, |side| map(side, |(&specie, &ν)| f64::powi(concentrations[specie], ν as i32)).product());
		let Rf = kf * ΠCνf;
		//dbg!(kr, ΠCνr);
		let Rr = kr * ΠCνr;
		//println!("{}", equation.format_with(" = ", |side, f| f(&side.format_with(" + ", |(&specie, ν), f| f(&format_args!("{} {}",ν,specie_names[specie]))))));
		//assert!(Rf >= Rr, Rf, Rr);
		let c = model.efficiency(T, &concentrations, kf);
		let R = Rf - Rr;
		let cR = c * R;
		for (specie, ν) in ν.enumerate() { dtω[specie] += ν * cR; }
	}
	let ref dtω = *dtω;

	/*{
		let (ref species, dtw) = {
			let X = std::ffi::CString::new(format!("H2:{}, O2:{}, AR:{}", concentrations[specie("H2").unwrap()]*W[specie("H2").unwrap()], concentrations[specie("O2").unwrap()]*W[specie("O2").unwrap()], C*W[specie("AR").unwrap()])).unwrap();
			let (mut len, mut species, mut dtw) = (0, std::ptr::null(), std::ptr::null());
			unsafe {
				cantera(/*rtol:*/ 1e-4, /*atol:*/ 1e-14, T, P, X.as_ptr(), &mut len, &mut species, &mut dtw);
				let species = std::slice::from_raw_parts(species, len);
				let dtw = std::slice::from_raw_parts(dtw, len);
				(box_collect(species.iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap())), dtw)
			}
		};
		print!("{}", species.zip(dtw.zip(species.map(|&s| specie(s).map(|k| dtω.get(k)).flatten()))).filter_map(|(k,(x,y))| y.map(|y| (k,x,y))).map(|(k,x,y)| format!("{} {:.1e} {:.1e}",k,x,y)).format(" "));
		std::io::Write::flush(&mut std::io::stdout()).unwrap();
		//std::process::abort();
	}*/

	let Cp = species.iter().map(|s| s.specific_heat_capacity(T));
	let rcp_ΣCCp = 1./concentrations.dot(Cp);
	let H = species.iter().map(|s| s.specific_enthalpy(T));
	let dtT = - rcp_ΣCCp * dtω.dot(H);
	let dtE = W.map(|w| 1.-w/W[S-1]).dot(dtω);
	let dtV = V * (dtT / T + T * ideal_gas_constant / P * dtE);
	let dtn = eval(dtω, |dtω| V*dtω);
	from_iter([dtT, dtV].chain(dtn))
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
		if i >= 2 && f64::abs(ρ - previous_ρ) <= 0.01*ρ.max(1./tmax) { dbg!(i); break; } // Early exit
		*yεv = eval!(y, dtεv; |y, dtεv| y + (ε / norm_dtεv) * dtεv);
	}
	dbg!(ρ);
	(yεv.sub(y), ρ * 1.2)
}

pub fn step(&self, rtol: f64, atol: f64, tmax: f64, P: f64, mut y: [f64; N]) -> [f64; N] {
	let max_steps = ((rtol / (10. * f64::EPSILON)).sqrt().round() as usize).max(2);
	let mut nstep = 0;
	let mut t = 0.;
	let ref mut dty = self.dt(P, &y);
	let (mut v, mut jacobian_spectral_radius) = self.power_iteration(P, tmax, &y, dty, dty);
	let mut dt = {
		let dt = (1./jacobian_spectral_radius).min(tmax);
		let ref dty1 = self.dt(P, &eval!(&y, &*dty; |y, dty| y + dt * dty));
		//(dt/(dt*error(iter::map!(&*dty, dty1, &y; |dty, dty1, y| (dty1 - dty) / (atol + rtol * y.abs())))) / 10.).min(tmax)
		(dt/(dt*error(zip!(&*dty, dty1, &y).map(|(dty, dty1, y):(&f64,&f64,&f64)| (dty1 - dty) / (atol + rtol * y.abs())))) / 10.).min(tmax)
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
		//let error = error(map!(&y,&y1,&*dty,dty1; |y,y1,dty,dty1| (0.8*(y1-y)+0.4*dt*(dty+dty1)/(atol + rtol*y.abs().max(y1.abs())))));
		let error = error(zip!(&y,&y1,&*dty,dty1).map(|(y,y1,dty,dty1):(&f64,&f64,&f64,&f64)| (0.8*(y1-y)+0.4*dt*(dty+dty1)/(atol + rtol*y.abs().max(y1.abs())))));
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
	::ron::de::from_str::<Map<Element, f64>>("#![enable(unwrap_newtypes)] {H: 1.008, O: 15.999, Ar: 39.95}").unwrap().map(|(e,g)| (e, g*1e-3/*kg/g*/)).collect()
});

impl<const S: usize, const S1: usize, const N: usize> Simulation<'t, S, S1, N> {
#[fehler::throws(anyhow::Error)] pub fn new(system: &'b [u8]) -> Self where 'b: 't {
	let ron::System{species: species_data, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;
	let Phase::IdealGas{species, state, ..} = phases.into_vec().into_iter().next().unwrap();
	let species : [_; S] = species.as_ref().try_into().unwrap();
	let molar_masses = eval(species, |s| species_data[s].composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum());
	let thermodynamics = eval(species, |s| { let Specie{thermodynamic:ron::NASA7{temperature_ranges,coefficients},..} = &species_data[s]; match temperature_ranges[..] {
		[_,Tsplit,_] if Tsplit == NASA7::T_split => NASA7(coefficients[..].try_into().unwrap()),
		[min, max] if min < NASA7::T_split && NASA7::T_split < max => NASA7([coefficients[0]; 2]),
		ref ranges => panic!("{:?}", ranges),
	}});
	let reactions = reactions.map(|self::ron::Reaction{ref equation, rate_constant, model}| {
		let atmospheric_pressure = 101325.;
		let equation = eval(equation, |e| box_collect(e.keys().map(|&key| species.iter().position(|&k| k==key).expect(key))).zip(box_collect(e.values().copied())));
		let specie_net_coefficients = generate(|specie| {
			let [reactant, product]:[_;2] = eval(&equation, |side| side.find(|(&s,_)| s==specie).map(|(_,&ν)| ν as i8).unwrap_or(0));
			(product-reactant) as f64
		});
		let PRν = f64::powf(atmospheric_pressure / ideal_gas_constant, specie_net_coefficients.iter().sum());
		Reaction{
			equation,
			rate_constant,
			model: {use self::ron::Model::*; match model {
				Elementary => Model::Elementary,
				ThreeBody{efficiencies} => Model::ThreeBody{efficiencies: eval(&species, |specie| *efficiencies.get(specie).unwrap_or(&1.))},
				Falloff{efficiencies, k0, troe} => Model::Falloff{efficiencies: eval(&species, |specie| *efficiencies.get(specie).unwrap_or(&1.)), k0, troe},
			}},
			specie_net_coefficients,
			//Σνf, Σνr,
			PRν
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

//impl<const S: usize> std::convert::From<State<S>> for Box<[Box<[f64]>]> { fn from(s: State<S>) -> Self { box [box [s.temperature] as Box<[_]>, box s.amounts] as Box<[_]> } }
