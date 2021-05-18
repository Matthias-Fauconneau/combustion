fn get_mut<T: Default, U>(cell: &std::cell::Cell<T>, f: impl FnOnce(&mut T) -> U) -> U {
	let mut value = cell.take();
	let result = f(&mut value);
	cell.set(value);
	result
}
fn implicit(u: &[f64]) -> Box<[f64]> { [u[0]].iter().chain(&u[2..]).copied().collect() } // one of P or V imply the other using ideal gas law
fn explicit(total_amount: f64, pressure_R: f64, u: &[f64]) -> Box<[f64]> { // Reconstructs P|V using ideal gas law
	let temperature = u[0];
	let volume = total_amount / pressure_R * temperature;
	[temperature, volume].iter().chain(&u[1..]).copied().collect()
}

use std::os::raw::c_char;
#[link(name = "cantera")]
extern "C" {
fn thermo_newFromFile(file_name: *const c_char, phase_name: *const c_char) -> i32;
fn thermo_nSpecies(n: i32) -> usize;
fn thermo_setTemperature(n: i32, t: f64) -> i32;
fn thermo_setMoleFractions(n: i32, len: usize, x: *const f64, norm: i32) -> i32;
fn thermo_getSpeciesName(n: i32, m: usize, len: usize, buffer: *mut c_char) -> i32;
fn thermo_setPressure(n: i32, p: f64) -> i32;
fn thermo_chemPotentials(n: i32, len: usize, murt: *mut f64);
fn kin_newFromFile(file_name: *const c_char, phase_name: *const c_char, reactingPhase: i32, neighbor0: i32, neighbor1: i32, neighbor2: i32, neighbor3: i32) -> i32;
fn kin_getEquilibriumConstants(n: i32, len: usize, kc: *mut f64) -> i32;
fn kin_getNetProductionRates(n: i32, len: usize, w_dot: *mut f64) -> i32;
fn kin_getFwdRatesOfProgress(n: i32, len: usize, rate: *mut f64) -> i32;
fn kin_getRevRatesOfProgress(n: i32, len: usize, rate: *mut f64) -> i32;
fn kin_getNetRatesOfProgress(n: i32, len: usize, rate: *mut f64) -> i32;
fn kin_getReactionString(n: i32, i: usize, len: usize, buffer: *mut c_char) -> i32;
}

use {fehler::throws, error::Error};
use itertools::Itertools;
use combustion::{*, reaction::{*, Property::*}};

#[throws] pub fn check(model: &model::Model, state: &State) {
	let (ref species_names, ref species) = Species::new(&model.species);
	let len = species.len();
	let reactions = iter::map(&model.reactions, |r| Reaction::new(species_names, r));
	let equations = iter::map(&reactions, |Reaction{reactants, products, model, ..}| {
		format!("{}", [reactants, products].iter().format_with(if let ReactionModel::Irreversible = model { " => " } else { " <=> " }, |side, f| {
			f(&side.iter().enumerate().filter(|(_,&ν)| ν > 0)
				.sorted_by_key(|&(k,_)| species_names[k])
				.format_with(" + ", |(specie, &ν), f| if ν > 1 { f(&format_args!("{} {}",ν,species_names[specie])) } else { f(&species_names[specie]) }))?;
			use ReactionModel::*; match model {
					Elementary|Irreversible => {},
					ThreeBody{..} => f(&" + M")?,
					PressureModification{..}|Falloff{..} => f(&" (+M)")?,
			};
			Ok(())
		}))
	});

	let file = std::ffi::CStr::from_bytes_with_nul(b"gri30.yaml\0").unwrap().as_ptr();
	let name = std::ffi::CStr::from_bytes_with_nul(b"gri30\0").unwrap().as_ptr();
	let phase = unsafe{thermo_newFromFile(file, name)};
	let cantera_species_names = iter::eval(species.len(), |k| {
		let mut specie = [0; 8];
		unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
		unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
	});
	assert!(unsafe{thermo_nSpecies(phase)} == len);
	assert!(state.amounts.len() == len && !state.amounts.iter().any(|&n| n<0.));
	let cantera_order = |o: &[f64]| (0..o.len()).map(|i| o[species_names.iter().position(|&s| s==cantera_species_names[i]).unwrap()]).collect::<Box<_>>();
	let kinetics = unsafe{kin_newFromFile(file, name, phase, 0, 0, 0, 0)};
	assert_eq!(equations, iter::eval(reactions.len(), |i| unsafe {
		let mut reaction = [0; 64]; kin_getReactionString(kinetics, i, reaction.len(), reaction.as_mut_ptr()); std::ffi::CStr::from_ptr(reaction.as_ptr()).to_str().unwrap().to_owned()
	}));
	let total_amount = state.amounts.iter().sum();
	let constant = state.constant();
	let state : StateVector<{Volume}> = state.into();
	let mut state = implicit(&state);
	let mut cvode = cvode::CVODE::new(&state);
	let mut time = 0.;
	let derivative = /*Derivative*/StateVector(std::iter::repeat(0.).take(2+len-1).collect());
	let derivative = std::cell::Cell::new(derivative);
	let (_, rate) = rate(species, &*reactions);
	// CVODE shim
	struct Derivative<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> {
		rate: Rate,
		constant: Constant<CONSTANT>,
		total_amount: f64,
		derivative: std::cell::Cell<StateVector::<CONSTANT>>
	}
	impl<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> FnOnce<(&[f64],)> for Derivative<Rate, CONSTANT> {
		type Output = Option<Box<[f64]>>;
		extern "rust-call" fn call_once(mut self, args: (&[f64],)) -> Self::Output { self.call_mut(args) }
	}
	impl<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> FnMut<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call_mut(&mut self, args: (&[f64],)) -> Self::Output { self.call(args) }
	}
	impl<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> Fn<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call(&self, (u,): (&[f64],)) -> Self::Output {
			let Self{rate, constant, total_amount, derivative} = self;
			get_mut(derivative, |mut derivative| {
				let pressure_R = constant.0 as f64;
				let mut debug = vec![f64::NAN; /*model.reactions.len()*/325*2].into_boxed_slice();
				//rate(*constant, &StateVector(map(&explicit(*total_amount, pressure, u), |&v| v.max(0.))), &mut derivative, &mut debug);
				rate(*constant, &StateVector(explicit(*total_amount, pressure_R, u)), &mut derivative, &mut debug);
				Some(implicit(&derivative.0))
			})
		}
	}
	let ref derivative = Derivative{
		rate,
		constant,
		total_amount,
		derivative
	};

	while std::hint::black_box(true) {
		let ref state_vector = StateVector(/*demote(&*/explicit(total_amount, constant.0 as f64, &state)/*)*/);
		let (ref equilibrium_constants, ref forward, ref reverse)/*(cantera_creation, cantera_destruction)*/ = {
			let state = State::new(total_amount, constant, &StateVector(iter::map(state_vector, |v| f64::max(0., *v))));
			unsafe{thermo_setMoleFractions(phase, state.amounts.len(), cantera_order(&state.amounts).as_ptr(), 1)}; // /!\ Needs to be set before pressure
			unsafe{thermo_setTemperature(phase, state.temperature)};
			unsafe{thermo_setPressure(phase, state.pressure_R * (K*NA))}; // /!\ Needs to be set after mole fractions
			let equilibrium_constants = {
				let mut equilibrium_constants = vec![0.; model.reactions.len()];
				unsafe{kin_getEquilibriumConstants(kinetics, model.reactions.len(), equilibrium_constants.as_mut_ptr())};
				equilibrium_constants.iter().zip(reactions.iter()).map(|(c,Reaction{Σnet, ..})| c*f64::powf(1e3, *Σnet as f64)).collect::<Box<_>>() // kmol -> mol
			};
			let forward = {
				let mut rates = vec![0.; model.reactions.len()];
				unsafe{kin_getFwdRatesOfProgress(kinetics, model.reactions.len(), rates.as_mut_ptr())};
				iter::map(&rates, |c| c*1000.) // kmol -> mol
			};
			let reverse = {
				let mut rates = vec![0.; model.reactions.len()];
				unsafe{kin_getRevRatesOfProgress(kinetics, model.reactions.len(), rates.as_mut_ptr())};
				iter::map(&rates, |c| c*1000.) // kmol -> mol
			};
			(equilibrium_constants, forward, reverse)
		};

		let mut debug = vec![f64::NAN; model.reactions.len()*2].into_boxed_slice();
		let _rate = {
			let ref rate = derivative.rate;
			let mut derivative = /*Derivative*/StateVector::<{Volume}>(vec![0.; 2+len-1].into_boxed_slice());
			rate(constant, state_vector, &mut derivative, &mut debug);
			derivative.0
		};
		let cR = &debug[0..model.reactions.len()];
		let rcp_equilibrium_constants = &debug[model.reactions.len()..2*model.reactions.len()];

		fn to_string(v: f64) -> String { if v == 0. { "0".to_owned() } else { format!("{:.0e}", v) } }

		if true {
			let T = state_vector[0];
			let a = iter::map(species.thermodynamics, |s| s.0[1]);
			let pressure_R = state_vector[1];
			let amounts = &state_vector[2..];
			fn dot(iter: impl IntoIterator<Item=(f64, f64)>) -> f64 { iter.into_iter().map(|(a,b)| a*b).sum() }
			use std::array::IntoIter;
			let [rcpT, T2, T3, T4] = [1./T, T*T, T*T*T, T*T*T*T];
			let exp_G_RT = iter::map(a[..len-1], |a| f64::exp(((a[0]-a[6]))+dot(
				IntoIter::new([(a[5], rcpT), (-a[0], f64::ln(T)), (-a[1]/2., T), ((1./3.-1./2.)*a[2], T2), ((1./4.-1./3.)*a[3], T3), ((1./5.-1./4.)*a[4], T4)]),
				)));
			for (((r, _e), (&rcpK, &cK)), (&_cR, (&forward, &reverse))) in reactions.iter().zip(equations.iter())
			.zip(rcp_equilibrium_constants.iter().zip(equilibrium_constants.iter()))
			.zip(cR.iter().zip(forward.iter().zip(reverse.iter()))) {
				fn product_of_exponentiations<T: Into<i16>>(iter: impl IntoIterator<Item=(T, f64)>) -> f64 { iter.into_iter().map(|(c,v)| f64::powi(v, c.into().into())).product() }
				let Reaction{reactants, products, rate_constant, net, Σnet, ..} = r;
				let P0_RT = NASA7::reference_pressure/T;
				let rcp_equilibrium_constant = product_of_exponentiations(net.iter().chain(&[-Σnet]).copied().zip(exp_G_RT.iter().chain(&[P0_RT]).copied()));
				fn fma(a: f64, b: f64, c: f64) -> f64 { a*b+c }
				fn arrhenius(RateConstant{preexponential_factor, temperature_exponent, activation_temperature}: RateConstant, T: f64) -> f64 {
					if [0.,-1.,1.,2.,4.,-2.].contains(&temperature_exponent) && activation_temperature == 0. {
						let A = preexponential_factor;
						if temperature_exponent == 0. { A }
						else if temperature_exponent == -1. { A/T }
						else if temperature_exponent == 1. { A*T }
						else if temperature_exponent == 2. { A*T*T }
						else if temperature_exponent == 4. { A*T*T*T*T }
						else if temperature_exponent == -2. { A/(T*T) }
						else { unreachable!() }
					} else {
						let logA = f64::ln(preexponential_factor);
						let βlogT𐊛logA = if temperature_exponent == 0. { logA } else { fma(temperature_exponent, f64::ln(T), logA) };
						let log_arrhenius = if activation_temperature == 0. { βlogT𐊛logA } else { fma(-activation_temperature, 1./T, βlogT𐊛logA) };
						f64::exp(log_arrhenius)
					}
				}
				let k_inf = arrhenius(*rate_constant, T);
				fn efficiency(model: &ReactionModel, T: f64, concentrations: &[f64], k_inf: f64) -> f64 {
					use ReactionModel::*; match model {
						Elementary|Irreversible => 1.,
						ThreeBody{efficiencies} => { dot(efficiencies.iter().copied().zip(concentrations.iter().copied())) },
						PressureModification{efficiencies, k0} => {
							let Pr = dot(efficiencies.iter().copied().zip(concentrations.iter().copied())) * arrhenius(*k0, T) / k_inf;
							Pr / (1.+ Pr)
						}
						Falloff{efficiencies, k0, troe} => {
							let Pr = dot(efficiencies.iter().copied().zip(concentrations.iter().copied())) * arrhenius(*k0, T) / k_inf;
							let model::Troe{A, T3, T1, T2} = *troe;
							fn rcp(x: f64) -> f64 { 1./x }
							let Fcent = fma(1.-A, f64::exp(-T*rcp(T3)), fma(A, f64::exp(-T*rcp(T1)), f64::exp(-1./T * T2)));
							let logFcent = f64::log2(Fcent);
							let c =fma(-0.67, logFcent, -0.4*f64::log2(10.));
							let N = fma(-1.27, logFcent, 0.75*f64::log2(10.));
							let logPr𐊛c = f64::log2(Pr)+c;
							let f1 = logPr𐊛c / fma(-0.14, logPr𐊛c, N);
							let F = f64::exp2(logFcent / fma(f1, f1, 1.));
							Pr / (1. + Pr) * F
						}
					}
				}

				let volume = constant.0;
				let rcpV = 1./volume;
				let total_concentration = pressure_R/T;
				let amounts = iter::map(amounts, |n| f64::max(0., *n));
				let concentrations = iter::map(&amounts, |&n| n*rcpV);
				let Ca = total_concentration - dot(std::iter::repeat(1.).zip(concentrations.iter().copied()));
				let ref concentrations = [&concentrations as &[_],&[Ca]].concat();
				let c = k_inf * efficiency(&r.model, T, concentrations, k_inf);
				let Rf = product_of_exponentiations(reactants.iter().copied().zip(concentrations.iter().copied()));
				{
					let (a, b) = (c*Rf, forward);
					use num::sign;
					assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
					assert!(num::relative_error(a, b) < 1e-3, "{:.3e}", num::relative_error(a, b));
				}
				if let ReactionModel::Irreversible = r.model {} else {
					let (a, b) = (1./rcpK, cK);
					use num::sign;
					assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
					assert!(num::relative_error(1./rcp_equilibrium_constant, 1./rcpK) < 1e-3, "{:e} {:e} {:e}", num::relative_error(1./rcpK, 1./rcp_equilibrium_constant), num::relative_error(1./rcpK, cK), num::relative_error(1./rcp_equilibrium_constant, cK));
					assert!(num::relative_error(a, b) < 1e-3, "{:.3e} {:e}", num::relative_error(a, b), num::relative_error(1./rcp_equilibrium_constant, cK));
				}
				{
					let Rr = if let ReactionModel::Irreversible = r.model { 0. } else {
						rcp_equilibrium_constant * product_of_exponentiations(products.iter().copied().zip(concentrations.iter().copied())) };
					let (a, b) = (c*Rr, reverse);
					use num::sign;
					assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
					assert!(num::relative_error(a, b) < 1e-5, "{:.3e}", num::relative_error(a, b));
				}
				let _ = rcp_equilibrium_constant;
			}
		}

		fn absolute_error(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(&a,&b)| f64::abs(a-b)).reduce(f64::max).unwrap() }
		fn relative_error(a: &[f64], b: &[f64]) -> f64 {
			a.iter().zip(b).map(|(&a,&b)|
				if a.abs() < 1e-3 && b.abs() < 1e-3 { 0. } else { num::relative_error(a,b) }
			).reduce(f64::max).unwrap()
		}

		let next_time = time + model.time_step;
		let mut steps = 0;
		while time < next_time {
			(time, state) = {
					let (time, state) = cvode.step(derivative, next_time, &state);
					(time, state.to_vec().into_boxed_slice())
			};
			steps += 1;
		}
		println!("{:.0} {:.0} {}", time*1e3, state[0], steps);
		assert_eq!(time, next_time);
	}
}
