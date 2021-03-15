//{std::ops::Deref, itertools::Itertools, super::*};
use super::*;

//pub fn promote(v: &[f32]) -> Box<[f64]> { v.iter().map(|&v| v as f64).collect() }

#[throws] pub fn check(model: Model, Simulation{/*species_names,*/ /*time_step,*/ /*mut*/ state, ..}: &Simulation) {
	let file = std::ffi::CStr::from_bytes_with_nul(b"gri30.yaml\0").unwrap().as_ptr();
	let name = std::ffi::CStr::from_bytes_with_nul(b"gri30\0").unwrap().as_ptr();
	let phase = unsafe{thermo_newFromFile(file, name)};
	assert!(model.len() == unsafe{thermo_nSpecies(phase)});
	/*let species = (0..model.len()).map(|k| {
		let mut specie = [0; 8];
		thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr());
		std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()
	});*/
	unsafe{thermo_setTemperature(phase, state.temperature)};
	unsafe{thermo_setPressure(phase, state.pressure)};
	assert!(state.amounts.len() == model.len());
	unsafe{thermo_setMoleFractions(phase, state.amounts.len(), state.amounts.as_ptr(), 1)};
	let kinetics = unsafe{kin_newFromFile(file, name, phase, 0, 0, 0, 0)};
	let mut net_production_rates = vec![0.; model.len()];
	unsafe{kin_getNetProductionRates(kinetics, model.len(), net_production_rates.as_mut_ptr())};
	println!("{}", net_production_rates);
	/*#[allow(unused_mut)] let mut time = 0.;
	const CONSTANT : Property = {use Property::*; Volume};
	//let mut cvode = cvode::CVODE::new(promote(((&state).into():StateVector<CONSTANT>).0.deref()).deref());
	while std::hint::black_box(true) {
		let next_time = time + *time_step;
		let (//equations, equilibrium_constants, [forward, reverse],
		ref cantera_rate,
		//ref cantera_state/*, _cantera_concentrations*/
		) /*: (_,_,_,Box<[f64]>,_)*/ = {
			//let initial_time = time;
			//let mole_proportions = state.amounts;
			let mole_proportions = format!("{}", species_names.iter().zip(state.amounts.deref()).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
			let mole_proportions = std::ffi::CString::new(mole_proportions)?;
			use std::ptr::null;
			let mut pressure = state.pressure / NA;
			//let volume = state.volume;
			let mut temperature = state.temperature;// / K;
			let (mut len,
						mut specie_names,
						mut net_productions_rates,
						//mut concentrations,
						//mut reactions_len, mut equations, mut equilibrium_constants, [mut forward, mut reverse]
						) = (
						0,
						null(),
						null(),
						//null(), 0, null(), null(), [null(); 2]
						);
			unsafe {
				cantera::reaction(&mut pressure, &mut temperature, mole_proportions.as_ptr(), &mut len, &mut specie_names,
																		//time-initial_time,
																		&mut net_productions_rates,
																		//&mut reactions_len, &mut equations, &mut equilibrium_constants, &mut forward, &mut reverse,
																		//next_time-initial_time, &mut concentrations
																		);
				let specie_names = iter::box_collect(std::slice::from_raw_parts(specie_names, len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
				let order = |o:&[_]| -> Box<[_]> { species_names.iter().map(|s| o[specie_names.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}", s, species_names))]).collect() };
				let net_productions_rates = order(std::slice::from_raw_parts(net_productions_rates, len)).iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>(); // kmol/m^3/s => mol/s [1m^3]
				//let concentrations = order(std::slice::from_raw_parts(concentrations, len)).iter().map(|c| c*1000.).collect(); // kmol/m^3 => mol/m^3
				(
					//iter::box_collect(std::slice::from_raw_parts(equations, reactions_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap())),
					//iter::box_collect(std::slice::from_raw_parts(equilibrium_constants, reactions_len).iter().zip(model.reactions.iter()).map(|(c,Reaction{Σnet, ..})| c*f64::powf(1e3, *Σnet as f64))),
					//[forward, reverse].map(|r| iter::box_collect(std::slice::from_raw_parts(r, reactions_len).iter().map(|c| c*1000.))),
					(net_productions_rates /* *volume*/),
					//State{temperature, pressure, volume, amounts: concentrations}
				)
			}
		};
		/*let cantera_reactions = iter::box_collect(iter::zip!(equations, equilibrium_constants, forward, reverse).into_iter());
		for (index, reaction) in cantera_reactions.iter().enumerate() {
			let net = reaction.2-reaction.3;
			if net != 0. { println!("{} {:.0e}", index, net); }
		}*/
		let (_, _, rate) = model.rate::<CONSTANT>();
		let mut derivative = /*Derivative*/StateVector::<CONSTANT>(std::iter::repeat(0.).take(model.len()).collect());
		let rate = {rate(state.constant(), &state.into(), &mut derivative); derivative};
		/*if false {
			let reactions = {
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
							let equation = format!("{}", [reactants, products].iter().format_with(" <=> ", |side, f| {
								//let other_species = ["Ar","H","H2","H2O","HO2","H2O2","O","O2","OH"];
								f(&side.iter().enumerate().filter(|(_,&ν)| ν > 0.)
									.sorted_by_key(|&(k,_)| species_names[k])
									.format_with(" + ", |(specie, &ν), f| if ν > 1. { f(&format_args!("{} {}",ν,species_names[specie])) } else { f(&species_names[specie]) }))?;
								use ReactionModel::*; match model {
										Elementary => {},
										ThreeBody{..} => f(&" + M")?,
										PressureModification{..}|Falloff{..} => f(&" (+M)")?,
								};
								Ok(())
							}));
							let log_kf = log_arrhenius(rate_constant, T);
							let c = model.efficiency(T, concentrations, log_kf);
							let mask = |mask, v| iter::zip!(mask, v).map(|(&mask, v):(_,&_)| if mask != 0. { *v } else { 0. });
							use iter::{into::IntoMap, vec::Dot};
							let Rf = f64::exp(reactants.dot(mask(reactants, log_concentrations)) + log_kf);
							let log_equilibrium_constant = -net.dot(G) + Σnet*logP0_RT;
							let Rr = f64::exp(products.dot(mask(products, log_concentrations)) + log_kf - log_equilibrium_constant);
							(equation, f64::exp(log_equilibrium_constant), c*Rf, c*Rr)
					}))
			};
			for (a, b) in reactions.iter().zip(cantera_reactions.iter()) {
					let ref e = a.0;
					//let b = other_reactions.into_iter().find(|(k,_,_,_)| k==&e).expect(&format!("{} {}", &e, other_reactions.iter().map(|b| b.0).format(" ")));
					assert_eq!(e, &b.0.replace(" =>"," <=>"));
					assert!(num::relative_error(a.1, b.1) < 0.002, "{}", num::relative_error(a.1, b.1));
					let a = a.2-a.3;
					let b = b.2-b.3;
					if num::relative_error(a,b) != 0. { println!("{:32} {:15.2e}", e, [a, b, num::relative_error(a,b)].iter().format(" ")); }
					use num::sign;
					assert!((sign(a)==sign(b) || (f64::max(f64::abs(a),f64::abs(b))<1e-17)) || (sign(a)==sign(b) && num::relative_error(a, b) < 0.));
			}
		}*/

		fn table(labels: &[&str], a: &[f64], b: &[f64]) -> Box<[([String; 4], usize)]> {
			labels.iter().zip(a.iter().zip(b)).filter(|(_,(&a,&b))| a != 0. || b != 0.).map(|(&header,(&a,&b))| {
				fn to_string(v: f64) -> String { if v == 0. { "0".to_owned() } else { format!("{:.0e}", v) } }
				let column = [header.to_owned(), to_string(a), to_string(b), to_string(num::relative_error(a,b))];
				let width = column.iter().map(|s| s.len()).max().unwrap();
				(column, width)
			}).collect()
		}
		fn print<const R: usize>(table: &[([String; R], usize)]) {
			for row in 0..R { println!("{}", table.iter().format_with(" ", |(c,width), f| f(&format_args!("{:width$}", c[row], width=width)))); }
		}

		let rate = promote(&rate[2..]);
		let rate = rate.deref();
		//let _ = &table(species_names.prefix(), rate, cantera_rate.prefix());
		let len = species_names.len();
		print(&table(&species_names[..len-1], rate, &cantera_rate[..len-1]));
		fn absolute_error(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(&a,&b)| f64::abs(a-b)).reduce(f64::max).unwrap() }
		fn relative_error(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(&a,&b)| num::relative_error(a,b)).reduce(f64::max).unwrap() }
		{
			let abs = absolute_error(rate, cantera_rate);
			let rel = relative_error(rate, cantera_rate);
			println!("rate {:e} {:e}", abs, rel);
			assert!(abs < 1e-8 && rel < 4e-4, "rate {:e} {:e}", abs, rel);
		}

		/*while time < next_time {
			let (next_time, next_state) = cvode.step(move |u| system.rate_and_jacobian::<CONSTANT>(state.constant::<CONSTANT>(), &State(*u)).map(|(rate, /*jacobian*/)| rate.0), next_time, &((&state).into():reaction::State<CONSTANT,S>)); //dbg!(time);
			(time, state) = (next_time, State::new(state.amounts.iter().sum(), state.constant::<CONSTANT>(), &reaction::State::<CONSTANT, S>(next_state)))
		}*/
		//let next_time = time;
		assert_eq!(time, next_time);
		println!("t {}", time);

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
	}*/
}
