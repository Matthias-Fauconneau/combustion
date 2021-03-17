pub fn map<T, U, F: Fn(&T)->U>(v: &[T], f: F) -> Box<[U]> { v.iter().map(f).collect() }
pub fn promote(v: &[f32]) -> Box<[f64]> { map(v, |&v| v as f64) }
pub fn demote(v: &[f64]) -> Box<[f32]> { map(v, |&v| v as f32) }
fn implicit(u: &[f64]) -> Box<[f64]> { [u[0]].iter().chain(&u[2..]).copied().collect() } // one of P or V imply the other using ideal gas law
fn explicit(total_amount: f64, pressure: f64, u: &[f64]) -> Box<[f64]> { // Reconstructs P|V using ideal gas law
	let temperature = u[0];
	let volume = K * total_amount / pressure * temperature;
	[temperature, volume].iter().chain(&u[1..]).copied().collect()
}

use super::{*, Property::*};
use itertools::Itertools;

#[throws] pub fn check(model: Model, Simulation{species_names, time_step, state, ..}: &Simulation) {
	let len = model.len();
	let equations = map(&model.reactions, |Reaction{reactants, products, model, ..}| {
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
	assert!(unsafe{thermo_nSpecies(phase)} == len);
	{
		let cantera_species_name = (0..len).map(|k| {
			let mut specie = [0; 8];
			unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
			unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
		}).collect::<Box<_>>();
		assert_eq!(&cantera_species_name.iter().map(String::as_str).collect::<Box<_>>(), species_names);
	}
	let kinetics = unsafe{kin_newFromFile(file, name, phase, 0, 0, 0, 0)};
	{
		let cantera_equations= (0..model.reactions.len()).map(|i| {
			let mut reaction = [0; 64];
			unsafe{kin_getReactionString(kinetics, i, reaction.len(), reaction.as_mut_ptr())};
			unsafe{std::ffi::CStr::from_ptr(reaction.as_ptr()).to_str().unwrap().to_owned()}
		}).collect::<Box<_>>();
		for (cantera, equation) in cantera_equations.iter().zip(equations.iter()) { assert_eq!(cantera, equation); }
		assert_eq!(cantera_equations, equations);
	}
	let total_amount = state.amounts.iter().sum();
	let volume = state.volume;
	let constant = state.constant();
	let mut state = implicit(&promote(&(state.into():StateVector<{Volume}>)));
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
	let (_, rate) = model.rate();
	/*let ref derivative = move |u| get_mut(derivative, |derivative| {
		let pressure = constant.0 as f64;
		rate(constant, &StateVector(map(&explicit(total_amount, pressure, u), |&v| (v as f32).max(0.))), &mut derivative);
		Some(implicit(&promote(&derivative.0)))
	});*/
	// CVODE shim
	struct Derivative<Rate: crate::Rate<CONSTANT>, const CONSTANT: Property> {
		rate: Rate,
		constant: Constant<CONSTANT>,
		total_amount: f64,
		derivative: std::cell::Cell<StateVector::<CONSTANT>>
	}
	impl<Rate: crate::Rate<CONSTANT>, const CONSTANT: Property> FnOnce<(&[f64],)> for Derivative<Rate, CONSTANT> {
		type Output = Option<Box<[f64]>>;
		extern "rust-call" fn call_once(mut self, args: (&[f64],)) -> Self::Output { self.call_mut(args) }
	}
	impl<Rate: crate::Rate<CONSTANT>, const CONSTANT: Property> FnMut<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call_mut(&mut self, args: (&[f64],)) -> Self::Output { self.call(args) }
	}
	impl<Rate: crate::Rate<CONSTANT>, const CONSTANT: Property> Fn<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call(&self, (u,): (&[f64],)) -> Self::Output {
			let Self{rate, constant, total_amount, derivative} = self;
			get_mut(derivative, |mut derivative| {
				let pressure = constant.0 as f64;
				let mut reactions = vec![0.; 325].into_boxed_slice();
				rate(*constant, &StateVector(map(&explicit(*total_amount, pressure, u), |&v| (v as f32).max(0.))), &mut derivative, &mut reactions);
				Some(implicit(&promote(&derivative.0)))
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
		let ref state_vector = StateVector(demote(&explicit(total_amount, constant.0 as f64, &state)));
		let (ref cantera_equilibrium_constants, ref _cantera_reactions, ref cantera_rates)/*(cantera_creation, cantera_destruction)*/ = {
			let state = State::new(total_amount, constant, state_vector);
			assert!(state.amounts.len() == len);
			unsafe{thermo_setMoleFractions(phase, state.amounts.len(), state.amounts.as_ptr(), 1)}; // /!\ Needs to be set before pressure
			unsafe{thermo_setTemperature(phase, state.temperature)};
			unsafe{thermo_setPressure(phase, state.pressure * NA)}; // /!\ Needs to be set after mole fractions
			let equilibrium_constants = {
				let mut equilibrium_constants = vec![0.; model.reactions.len()];
				unsafe{kin_getEquilibriumConstants(kinetics, model.reactions.len(), equilibrium_constants.as_mut_ptr())};
				equilibrium_constants.iter().zip(model.reactions.iter()).map(|(c,Reaction{Σnet, ..})| c*f64::powf(1e3, *Σnet as f64)).collect::<Box<_>>() // kmol -> mol
			};
			let reaction_rates = {
				let mut rates = vec![0.; model.reactions.len()];
				unsafe{kin_getNetRatesOfProgress(kinetics, model.reactions.len(), rates.as_mut_ptr())};
				rates.iter().map(|c| c*1000.).collect::<Box<_>>() // kmol -> mol
			};
			let rates = {
				let mut rates = vec![0.; len];
				unsafe{kin_getNetProductionRates(kinetics, len, rates.as_mut_ptr())};
				rates.iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>() // kmol -> mol
			};
			/*let creation = {
				let mut rates = vec![0.; len];
				unsafe{kin_getCreationRates(kinetics, len, rates.as_mut_ptr())};
				rates.iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>() // kmol -> mol
			}
			let destruction = {
				let mut rates = vec![0.; len];
				unsafe{kin_getDestructionRates(kinetics, len, destruction_rates.as_mut_ptr())};
				rates.iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>() // kmol -> mol
			}
			(creation, destruction)*/
			(equilibrium_constants, reaction_rates, rates)
		};

		let mut rcp_equilibrium_constants = vec![f32::NAN; model.reactions.len()].into_boxed_slice();
		//let mut reactions_rates = vec![0.; model.reactions.len()].into_boxed_slice();
		//let mut [creation, destruction] = [vec![0.; len].into_boxed_slice(); 2];
		let rate = {
			let ref rate = derivative.rate;
			let mut derivative = /*Derivative*/StateVector::<{Volume}>(vec![0.; 2+len-1].into_boxed_slice());
			rate(constant, state_vector, &mut derivative, &mut rcp_equilibrium_constants);
			//rate(constant, state_vector, &mut derivative, &mut reactions_rates);
			//rate(constant, state_vector, &mut derivative, [&mut creation, &mut destruction]);
			derivative.0
		};

		fn to_string(v: f64) -> String { if v == 0. { "0".to_owned() } else { format!("{:.0e}", v) } }

		if false {
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
			let Model{species: Species{thermodynamics, ..}, reactions, ..} = &model;
			let a = thermodynamics.iter().map(|s| s.0[1]).collect(): Box<_>;
			let T = state_vector[0];
			let [rcpT, logT, T2, T3, T4] = [1./T, f32::log2(T), T*T, T*T*T, T*T*T*T];
			fn dot(iter: impl IntoIterator<Item=(f64, f32)>) -> f32 { iter.into_iter().map(|(a,b)| (a as f32)*b).sum() }
			use std::array::IntoIter;
			use std::f64::consts::LN_2;
			let exp_G_RT = a[..len-1].iter().map(|a| f32::exp2(((a[0]-a[6])/LN_2) as f32+dot(
				IntoIter::new([(a[5]/LN_2, rcpT), (-a[0], logT), (-a[1]/2./LN_2, T), ((1./3.-1./2.)*a[2]/LN_2, T2), ((1./4.-1./3.)*a[3]/LN_2, T3), ((1./5.-1./4.)*a[4]/LN_2, T4)]),
				))).collect():Box<_>;
			for ((r, e), (&a, &b)) in reactions.iter().zip(equations.iter()).zip(rcp_equilibrium_constants.iter().zip(cantera_equilibrium_constants.iter())) {
				if let ReactionModel::Irreversible = r.model { continue; }
				let a = 1./a as f64;
				if num::relative_error(a,b) != 0. { println!("{:32} {:15.2e}", e, [a, b, num::relative_error(a,b)].iter().format(" ")); }
				use num::sign;
				assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
				//fn product_of_exponentiations(iter: impl IntoIterator<Item=(i8, f32)>) -> f32 { iter.into_iter().map(|(c,v)| f32::powi(v, c.into())).product() }
				fn product_of_exponentiations(iter: impl IntoIterator<Item=(i8, f32)>) -> f64 { iter.into_iter().map(|(c,v)| f64::powi(v as f64, c.into())).product() }
				let Reaction{net, ..} = r;
				let rcp_equilibrium_constant = product_of_exponentiations(net.iter().copied().zip(exp_G_RT.iter().copied()));
				//let rcp_equilibrium_constant = product_of_exponentiations(net.iter().chain(&[-Σnet]).copied().zip(exp_G_RT.iter().chain(&[P0_RT]).copied()), None, C, f).unwrap();
				assert!(num::relative_error(a, b) < 0.03, "{:?} {:e} {:e}", net.iter().zip(exp_G_RT.iter()).filter(|(&c,_)| c != 0).map(|(c,&v)| (c, to_string(v as f64))).format(" "), rcp_equilibrium_constant, 1./rcp_equilibrium_constant);
			}
			/*for (e, (&a, &b)) in reactions.iter().zip(reactions_rates.iter().zip(cantera_reactions.iter())) {
				let a = a as f64;
				if num::relative_error(a,b) != 0. { println!("{:32} {:15.2e}", e, [a, b, num::relative_error(a,b)].iter().format(" ")); }
				use num::sign;
				assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
				assert!(num::relative_error(a, b) < 0.002, "{:.1e} {:.1e} {}", a, b, num::relative_error(a, b));
			}*/
		}

		let rate = &rate[2..];
		let rate = promote(&rate);
		let ref rate = rate.iter().map(|dtn| dtn/volume).collect::<Box<_>>();
		//let ref rate = rate.iter().map(|&dtn| if dtn.abs() < 1e-29 { 0. } else { dtn }).collect::<Box<_>>(); // Correct max relative error
		fn absolute_error(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(&a,&b)| f64::abs(a-b)).reduce(f64::max).unwrap() }
		fn relative_error(a: &[f64], b: &[f64]) -> f64 {
			a.iter().zip(b).map(|(&a,&b)|
				if a.abs() < 1e-2 && b.abs() < 1e-2 { 0. } else { num::relative_error(a,b) }
			).reduce(f64::max).unwrap()
		}
		let abs = absolute_error(rate, cantera_rates);
		let rel = relative_error(rate, cantera_rates);
		println!("{} {:.5} {:e} {:e}", time*1e3, state[0], abs, rel);
		//last_time = time;
		if !(abs < 7e-3 && rel < 1e-2) {
			fn table(labels: &[&str], a: &[f64], b: &[f64]) -> Box<[([String; 5], usize)]> {
				labels.iter().zip(a.iter().zip(b)).filter(|(_,(&a,&b))| a.abs() > 1e-29 || b.abs() > 1e-29).map(|(&header,(&a,&b))| {
					let column = [header.to_owned(), to_string(a), to_string(b), to_string(num::abs(a-b)), to_string(num::relative_error(a,b))];
					let width = column.iter().map(|s| s.len()).max().unwrap();
					(column, width)
				}).collect()
			}
			fn print<const R: usize>(table: &[([String; R], usize)]) {
				for row in 0..R { println!("{}", table.iter().format_with(" ", |(c,width), f| f(&format_args!("{:width$}", c[row], width=width)))); }
			}
			print(&table(&species_names[..len-1], rate, &cantera_rates[..len-1]));
		}
		assert!(abs < 7e-3 && rel < 1e-2, "{:e} {:e}", abs, rel);

		let next_time = time + time_step;
		let mut steps = 0;
		while time < next_time {
			(time, state) = {
					let (time, state) = cvode.step(derivative, next_time, &state);
					(time, state.to_vec().into_boxed_slice())
			};
			steps += 1;
		}
		println!("{}", steps);
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
