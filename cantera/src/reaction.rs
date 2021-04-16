#![allow(dead_code)]
pub fn map<T, U, F: Fn(&T)->U>(v: &[T], f: F) -> Box<[U]> { v.iter().map(f).collect() }
//pub fn promote(v: &[f32]) -> Box<[f64]> { map(v, |&v| v as f64) }
//pub fn demote(v: &[f64]) -> Box<[f32]> { map(v, |&v| v as f32) }
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
	let reactions = map(&model.reactions, |r| Reaction::new(species_names, r));
	let equations = map(&reactions, |Reaction{reactants, products, model, ..}| {
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
	{
		let cantera_equations= (0..reactions.len()).map(|i| {
			let mut reaction = [0; 64];
			unsafe{kin_getReactionString(kinetics, i, reaction.len(), reaction.as_mut_ptr())};
			unsafe{std::ffi::CStr::from_ptr(reaction.as_ptr()).to_str().unwrap().to_owned()}
		}).collect::<Box<_>>();
		for (cantera, equation) in cantera_equations.iter().zip(equations.iter()) { assert_eq!(cantera, equation); }
		assert_eq!(cantera_equations, equations);
	}
	let total_amount = state.amounts.iter().sum();
	//let volume = state.volume;
	let constant = state.constant();
	let mut state = implicit(/*&promote(*/&(state.into():StateVector<{Volume}>)/*)*/);
	let mut cvode = cvode::CVODE::new(&state);
	let mut time = 0.;
	let derivative = /*Derivative*/StateVector(std::iter::repeat(0.).take(2+len-1).collect());
	let derivative = std::cell::Cell::new(derivative);
	fn get_mut<T: Default, U>(cell: &std::cell::Cell<T>, f: impl FnOnce(&mut T) -> U) -> U {
		let mut value = cell.take();
		let result = f(&mut value);
		cell.set(value);
		result
	}
	let (_, rate) = rate(species, &*reactions);
	/*let ref derivative = move |u| get_mut(derivative, |derivative| {
		let pressure = constant.0 as f64;
		rate(constant, &StateVector(map(&explicit(total_amount, pressure, u), |&v| (v as f32).max(0.))), &mut vec![f64::NAN; model.reactions.len()].into_boxed_slice());
		Some(implicit(&promote(&derivative.0)))
	});*/
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

	//let mut last_time = time;
	while std::hint::black_box(true) {
		let ref state_vector = StateVector(/*demote(&*/explicit(total_amount, constant.0 as f64, &state)/*)*/);
		let (ref equilibrium_constants, ref forward, ref reverse)/*(cantera_creation, cantera_destruction)*/ = {
			let state = State::new(total_amount, constant, &StateVector(map(state_vector, |v| f64::max(0., *v))));
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
				map(&rates, |c| c*1000.) // kmol -> mol
			};
			let reverse = {
				let mut rates = vec![0.; model.reactions.len()];
				unsafe{kin_getRevRatesOfProgress(kinetics, model.reactions.len(), rates.as_mut_ptr())};
				map(&rates, |c| c*1000.) // kmol -> mol
			};
			//let order = |o: &[f64]| (0..species.len()).map(|i| o[cantera_species_names.iter().position(|s| s==species_names[i]).unwrap()]).collect::<Box<_>>();
			/*let rates = {
				let mut rates = vec![0.; len];
				unsafe{kin_getNetProductionRates(kinetics, len, rates.as_mut_ptr())};
				order rates.iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>() // kmol -> mol
			};*/
			/*let creation = {
				let mut rates = vec![0.; len];
				unsafe{kin_getCreationRates(kinetics, len, rates.as_mut_ptr())};
				order rates.iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>() // kmol -> mol
			}
			let destruction = {
				let mut rates = vec![0.; len];
				unsafe{kin_getDestructionRates(kinetics, len, destruction_rates.as_mut_ptr())};
				order rates.iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>() // kmol -> mol
			}
			(creation, destruction)*/
			(equilibrium_constants, forward, reverse)
		};

		let mut debug = vec![f64::NAN; model.reactions.len()*2].into_boxed_slice();
		//let mut reactions_rates = vec![0.; model.reactions.len()].into_boxed_slice();
		//let mut [creation, destruction] = [vec![0.; len].into_boxed_slice(); 2];
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
			/*let reactions = {
					let Model{species: Species{thermodynamics, ..}, reactions, ..} = &model;
					let State{volume, temperature, amounts, ..} = state;
					let T = temperature;
					use num::log;
					let logP0_RT = log(NASA7::reference_pressure) - log(T);
					use iter::{vec::eval, eval};
					let ref H = eval(thermodynamics, |s| s.specific_enthalpy(T));
					let ref H_T = eval(H.prefix(), |H| H/T);
					let ref G = eval!(thermodynamics.prefix(), H_T; |s, h_T| h_T - s.specific_entropy(T)); // (H-TS)/RT
					let ref concentrations = eval(amounts, |n| n / volume);
					let ref log_concentrations = eval(concentrations, |&c| f64::ln(c));
					iter::box_collect(reactions.iter().map(move |Reaction{reactants, products, rate_constant, model, net, Σnet, ..}| {
							let log_kf = log_arrhenius(rate_constant, T);
							let c = model.efficiency(T, concentrations, log_kf);
							let mask = |mask, v| iter::zip!(mask, v).map(|(&mask, v):(_,&_)| if mask != 0. { *v } else { 0. });
							use iter::{into::IntoMap, vec::Dot};
							let Rf = f64::exp(reactants.dot(mask(reactants, log_concentrations)) + log_kf);
							let log_equilibrium_constant = -net.dot(G) + Σnet*logP0_RT;
							let Rr = f64::exp(products.dot(mask(products, log_concentrations)) + log_kf - log_equilibrium_constant);
							(equation, f64::exp(log_equilibrium_constant), c*Rf, c*Rr)
					}))
			};*/
			let T = state_vector[0];
			//assert!(T >= NASA7::T_split, "{} {}", T, NASA7::T_split);
			let a = species.thermodynamics.iter().map(|s| s.0[1]).collect(): Box<_>;
			let pressure_R = state_vector[1];
			//assert_eq!(pressure_R, 101325./(K*NA));
			let amounts = &state_vector[2..];
			fn dot(iter: impl IntoIterator<Item=(f64, f64)>) -> f64 { iter.into_iter().map(|(a,b)| a*b).sum() }
			use std::array::IntoIter;
			let [rcpT, T2, T3, T4] = [1./T, T*T, T*T*T, T*T*T*T];
			/*let G0_RT = a[..len-1].iter().map(|a| ((a[0]-a[6]))+dot(
				IntoIter::new([(a[5], rcpT), (-a[0], f64::ln(T)), (-a[1]/2., T), ((1./3.-1./2.)*a[2], T2), ((1./4.-1./3.)*a[3], T3), ((1./5.-1./4.)*a[4], T4)]),
				)).collect():Box<_>;
			let RT = NA*K*T;
			let ref standard_chemical_potentials = map(&G0_RT, |g| g*RT+f64::ln(pressure_R / NASA7::reference_pressure)*RT);
			dbg!(standard_chemical_potentials);
			//let ref cantera_standard_chemical_potentials = RT() * log(xx);
			//dbg!(chemical_potentials);
			//for (&a,&b) in cantera_chemical_potentials.iter().zip(chemical_potentials.iter()) { assert!(num::relative_error(a,b) < 1e-100, "{:e}", num::relative_error(a,b)); }*/
			let exp_G_RT = a[..len-1].iter().map(|a| f64::exp(((a[0]-a[6]))+dot(
				IntoIter::new([(a[5], rcpT), (-a[0], f64::ln(T)), (-a[1]/2., T), ((1./3.-1./2.)*a[2], T2), ((1./4.-1./3.)*a[3], T3), ((1./5.-1./4.)*a[4], T4)]),
				))).collect():Box<_>;
			for (((r, e), (&rcpK, &cK)), (&_cR, (&forward, &reverse))) in reactions.iter().zip(equations.iter())
			.zip(rcp_equilibrium_constants.iter().zip(equilibrium_constants.iter()))
			.zip(cR.iter().zip(forward.iter().zip(reverse.iter()))) {
				//if let ReactionModel::Irreversible = r.model { continue; }
				//if num::relative_error(1./rcpK,cK) != 0. { println!("{:32} {:15.2e}", e, [1./rcpK, cK, num::relative_error(1./rcpK,cK)].iter().format(" ")); }
				//use num::sign;
				//assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
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
				let amounts = map(amounts, |n| f64::max(0., *n));
				let concentrations = map(&amounts, |&n| n*rcpV);
				let Ca = total_concentration - dot(std::iter::repeat(1.).zip(concentrations.iter().copied()));
				let ref concentrations = [&concentrations as &[_],&[Ca]].concat();
				let c = k_inf * efficiency(&r.model, T, concentrations, k_inf);
				let Rf = product_of_exponentiations(reactants.iter().copied().zip(concentrations.iter().copied()));
				{
					let (a, b) = (c*Rf, forward);
					/*if num::relative_error(a,b) != 0.*/ { println!("{:32} fwd {:15.2e}", e, [a, b, num::relative_error(a,b)].iter().format(" ")); }
					use num::sign;
					assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
					assert!(num::relative_error(a, b) < 1e-3, "{:.3e}", num::relative_error(a, b));
					/*assert!(num::relative_error(a, b) < 1e-3, "{:?} {:.3e} {:.3e} {:.3e} {:.3e}", net.iter().zip(exp_G_RT.iter()).filter(|(&c,_)| c != 0).map(|(c,&v)| (c, to_string(v as f64))).format(" "),
					1./rcp_equilibrium_constant, b, num::relative_error(b, 1./rcp_equilibrium_constant), num::relative_error(cR, ccR));*/
				}
				if let ReactionModel::Irreversible = r.model {} else {
					let (a, b) = (1./rcpK, cK);
					if num::relative_error(a,b) != 0. { println!("{:32} K {:15.3e}", e, [a, b, num::relative_error(a,b)].iter().format(" ")); }
					use num::sign;
					assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
					assert!(num::relative_error(1./rcp_equilibrium_constant, 1./rcpK) < 1e-3, "{:e} {:e} {:e}", num::relative_error(1./rcpK, 1./rcp_equilibrium_constant), num::relative_error(1./rcpK, cK), num::relative_error(1./rcp_equilibrium_constant, cK));
					assert!(num::relative_error(a, b) < 1e-3, "{:.3e} {:e}", num::relative_error(a, b), num::relative_error(1./rcp_equilibrium_constant, cK));
					/*assert!(num::relative_error(a, b) < 1e-3, "{:?} {:.3e} {:.3e} {:.3e} {:.3e}", net.iter().zip(exp_G_RT.iter()).filter(|(&c,_)| c != 0).map(|(c,&v)| (c, to_string(v as f64))).format(" "),
					1./rcp_equilibrium_constant, b, num::relative_error(b, 1./rcp_equilibrium_constant), num::relative_error(cR, ccR));*/
				}
				{
					let Rr = if let ReactionModel::Irreversible = r.model { 0. } else {
						rcp_equilibrium_constant * product_of_exponentiations(products.iter().copied().zip(concentrations.iter().copied())) };
					let (a, b) = (c*Rr, reverse);
					if num::relative_error(a,b) != 0. { println!("{:32} rev {:15.2e}", e, [a, b, num::relative_error(a,b)].iter().format(" ")); }
					use num::sign;
					assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
					assert!(num::relative_error(a, b) < 1e-5, "{:.3e}", num::relative_error(a, b));
					/*assert!(num::relative_error(a, b) < 1e-3, "{:?} {:.3e} {:.3e} {:.3e} {:.3e}", net.iter().zip(exp_G_RT.iter()).filter(|(&c,_)| c != 0).map(|(c,&v)| (c, to_string(v as f64))).format(" "),
					1./rcp_equilibrium_constant, b, num::relative_error(b, 1./rcp_equilibrium_constant), num::relative_error(cR, ccR));*/
				}
				let _ = rcp_equilibrium_constant;
				//assert!(num::relative_error(a, b) < 1e-3, "{:?} {:.3e} {:.3e}", net.iter().zip(exp_G_RT.iter()).filter(|(&c,_)| c != 0).map(|(c,&v)| (c, to_string(v as f64))).format(" "), num::relative_error(b, 1./rcp_equilibrium_constant), num::relative_error(cR, ccR));
			}
			/*println!("cR");
			for ((_r, e), (&a, &b)) in model.reactions.iter().zip(equations.iter()).zip(cR.iter().zip(_cantera_reactions.iter())) {
				if num::relative_error(a,b) != 0. { println!("{:32} {:15.2e}", e, [a, b, num::relative_error(a,b)].iter().format(" ")); }
				use num::sign;
				assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
				assert!(num::relative_error(a, b) < 1e-3);
			}*/
			/*for (e, (&a, &b)) in reactions.iter().zip(reactions_rates.iter().zip(cantera_reactions.iter())) {
				let a = a as f64;
				if num::relative_error(a,b) != 0. { println!("{:32} {:15.2e}", e, [a, b, num::relative_error(a,b)].iter().format(" ")); }
				use num::sign;
				assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
				assert!(num::relative_error(a, b) < 0.002, "{:.1e} {:.1e} {}", a, b, num::relative_error(a, b));
			}*/
		}

		//let rate = &rate[2..];
		//let rate = promote(&rate);
		//let ref rate = rate.iter().map(|dtn| dtn/volume).collect::<Box<_>>();
		fn absolute_error(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(&a,&b)| f64::abs(a-b)).reduce(f64::max).unwrap() }
		fn relative_error(a: &[f64], b: &[f64]) -> f64 {
			a.iter().zip(b).map(|(&a,&b)|
				if a.abs() < 1e-3 && b.abs() < 1e-3 { 0. } else { num::relative_error(a,b) }
			).reduce(f64::max).unwrap()
		}
		/*let abs = absolute_error(rate, cantera_rates);
		let rel = relative_error(rate, cantera_rates);
		//println!("{} {:.5} {:e} {:e}", time*1e3, state[0], abs, rel);
		//last_time = time;
		if !(abs < 7e-5 && rel < 1e-3) {
			/*fn table(labels: &[&str], a: &[f64], b: &[f64]) -> Box<[([String; 5], usize)]> {
				labels.iter().zip(a.iter().zip(b)).filter(|(_,(&a,&b))| a.abs() > 1e-29 || b.abs() > 1e-29).map(|(&header,(&a,&b))| {
					let column = [header.to_owned(), to_string(a), to_string(b), to_string(num::abs(a-b)), to_string(num::relative_error(a,b))];
					let width = column.iter().map(|s| s.len()).max().unwrap();
					(column, width)
				}).collect()
			}
			fn print<const R: usize>(table: &[([String; R], usize)]) {
				for row in 0..R { println!("{}", table.iter().format_with(" ", |(c,width), f| f(&format_args!("{:width$}", c[row], width=width)))); }
			}
			print(&table(&species_names[..len-1], rate, &cantera_rates[..len-1]));*/
			for (specie, (&rate, &cantera)) in species_names[..len-1].iter().zip(rate.iter().zip(cantera_rates[..len-1].iter())).filter(|(_,(&a,&b))| a.abs() > 1e-29 || b.abs() > 1e-29) {
				println!("{:6} {:5.0e} {:5.0e} {:5.0e} {:5.0e}", specie, rate, cantera, num::abs(rate-cantera), num::relative_error(rate, cantera));
			}
		}
		assert!(abs < 7e-3 && rel < 1e-3, "{:e} {:e}", abs, rel);*/

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
		//println!("{:.0} {:.0} {:.0e} {:.0e} {}", time*1e3, state[0], abs, rel, steps);
		assert_eq!(time, next_time);
		/*println!("t {}", time);
		(time, state) = {
			let (time, state) = cvode.step(derivative, time+time_step, &state);
			(time, state.to_vec().into_boxed_slice())
		}*/

		/*{
			println!("T {} {} {:e}", state.temperature, cantera_state.temperature, num::relative_error(state.temperature, cantera_state.temperature));
			print(&table(species_names, &state.amounts, &cantera_state.amounts));
			{
				let abs = absolute_error(&state.amounts, &cantera_state.amounts);
				let rel = relative_error(&state.amounts, &cantera_state.amounts);
				println!("state {:e} {:e}", abs, rel);
				assert!(abs < 1e-8 || rel < 0., "state {:e} {:e}", abs, rel);
			}
		}*/

		//state = *cantera_state; // Check rates along cantera trajectory
		//break;
	}
}
