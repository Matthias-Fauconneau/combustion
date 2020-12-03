//#![allow(incomplete_features)]#![feature(const_generics, const_evaluatable_checked,
#![feature(type_ascription, array_map, non_ascii_idents)]#![allow(mixed_script_confusables,non_snake_case)]
extern "C" {
fn cantera(relative_tolerance: f64, absolute_tolerance: f64, temperature: &mut f64, pressure: &mut f64, mole_proportions: *const std::os::raw::c_char, time_step: f64,
									species_len: &mut usize, species: &mut *const *const std::os::raw::c_char, net_production_rates: &mut *const f64, concentrations: &mut *const f64,
									reactions_len: &mut usize, equations: &mut *const *const std::os::raw::c_char, equilibrium_constants: &mut *const f64, forward: &mut *const f64, reverse: &mut *const f64);
}

#[fehler::throws(Box<dyn std::error::Error>)] fn main() {
	use combustion::*;
	let system = std::fs::read("H2+O2.ron")?;
	pub const S : usize = 9; // Number of species
	let Simulation{species, system, state: State{temperature, ref amounts}, volume, pressure, time_step, ..} = Simulation/*::<S>*/::new(&system)?;
	//type Vec = Simulation::<S>::Vec;
	let mut state /*: [_; Simulation::<S>::state_vector_len]*/ = {
			use {std::convert::TryInto, iter::{into::IntoChain, array_from_iter as from_iter}};
			from_iter([temperature,volume].chain(amounts[..S-1].try_into().unwrap():[_;S-1]))
	};
	/*for _ in 0..2 {
	use itertools::Itertools;
	let (equations, [equilibrium_constants, forward, reverse], ref other_net_production_rates, ref other_concentrations) = {
		let (T, _, ref concentrations) /*: (_,_,[_;S])*/ = System::state(pressure, &state);
		let mole_proportions = format!("{}", species.iter().zip(concentrations).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
		let mole_proportions = std::ffi::CString::new(mole_proportions).unwrap();
		use std::ptr::null;
		let (mut T, mut pressure, mut species_len, mut specie_names, mut net_production_rates, mut concentrations, mut reactions_len, mut equations, [mut equilibrium_constants, mut forward, mut reverse]) =
					(T, pressure, 0, null(), null(), null(), 0, null(), [null(); 3]);
		unsafe {
			cantera(/*relative_tolerance:*/ 1e-8, /*absolute_tolerance:*/ 1e-14, &mut T, &mut pressure, mole_proportions.as_ptr(), time_step,
				&mut species_len, &mut specie_names,  &mut net_production_rates, &mut concentrations,
				&mut reactions_len, &mut equations, &mut equilibrium_constants, &mut forward, &mut reverse);
			let equations = iter::box_collect(std::slice::from_raw_parts(equations, reactions_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
			let equilibrium_constants = std::slice::from_raw_parts(equilibrium_constants, reactions_len);
			let equilibrium_constants = iter::box_collect(equilibrium_constants.iter().zip(system.reactions.iter()).map(|(c,Reaction{Σnet, ..})| c*f64::powf(1e3, *Σnet)));
			let [forward, reverse] = [forward, reverse].map(|r| iter::box_collect(std::slice::from_raw_parts(r, reactions_len).iter().map(|c| c*1000.)));
			let specie_names = iter::box_collect(std::slice::from_raw_parts(specie_names, species_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
			let order = |o:Box<[_]>| iter::vec::eval(species, |s| o[specie_names.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}",s,species))]);
			let net_production_rates = iter::box_collect(std::slice::from_raw_parts(net_production_rates, species_len).iter().map(|c| c*1000.)); // kmol/m^3/s => mol/s [1m^3]
			let concentrations = iter::box_collect(std::slice::from_raw_parts(concentrations, species_len).iter().map(|c| c*1000.)); // kmol/m^3 => mol/m^3
			(equations, [equilibrium_constants, forward, reverse], std::convert::TryInto::try_into(&order(net_production_rates)[..S-1]).unwrap(), order(concentrations))
		}
	};
	if true {
	let other_reactions = iter::box_collect(iter::zip!(equations, equilibrium_constants, forward, reverse).into_iter());
	let reactions = {
		let specie_names = species;
		let System{thermodynamics: species, reactions, ..} = &system;
		let (T, _, ref concentrations) : (_,_,[_;S]) = System::state(pressure, &state);
		println!("{}", T);
		let logP0_RT = f64::ln(NASA7::reference_pressure/ideal_gas_constant) - f64::ln(T);
		let ref H_T = iter::vec::eval(species, |s| s.dimensionless_specific_enthalpy_T(T));
		let ref G = iter::eval!(species, H_T; |s, h_T| h_T - s.dimensionless_specific_entropy(T)); // (H-TS)/RT
		let log_concentrations = iter::vec::eval(concentrations, |&c| f64::ln(c));
		iter::box_collect(reactions.iter().map(move |Reaction{reactants, products, rate_constant, model, net, Σnet, ..}| {
			let equation = format!("{}", [reactants, products].iter().format_with(" <=> ", |side, f| {
				let other_species = ["Ar","H","H2","H2O","HO2","H2O2","O","O2","OH"];
				f(&
					//(&side).into_iter()
					Iterator::zip(side.a.into_iter(), side.b.into_iter())
					.sorted_by_key(|(&k,_)| other_species.iter().position(|s| s==&specie_names[k]).unwrap()).format_with(" + ", |(&specie, &ν), f| if ν > 1. { f(&format_args!("{} {}",ν,specie_names[specie].to_uppercase())) } else { f(&specie_names[specie].to_uppercase()) }))?;
				use Model::*; match model {
					Elementary => {},
					ThreeBody{..} => f(&" + M")?,
					Falloff{..} => f(&" (+M)")?,
				};
				Ok(())
			}));
			let log_kf = log_arrhenius(rate_constant, T);
			let Rf = f64::exp(reactants.dot(log_concentrations) + log_kf);
			let log_equilibrium_constant = -net.dot(G) + Σnet*logP0_RT;
			let Rr = f64::exp(products.dot(log_concentrations) + log_kf - log_equilibrium_constant);
			let [Rf,Rr] = [Rf,Rr].map(|R| model.efficiency(T, &concentrations, log_kf) * R);
			(equation, f64::exp(log_equilibrium_constant), Rf, Rr)
		}))
	};
	for (a, b) in reactions.iter().zip(other_reactions.iter()) {
		let ref e = a.0;
		//let b = other_reactions.into_iter().find(|(k,_,_,_)| k==&e).expect(&format!("{} {}", &e, other_reactions.iter().map(|b| b.0).format(" ")));
		assert!(e == b.0);
		let a = a.2-a.3;
		let b = b.2-b.3;
		println!("{:30} {:15.2e}", e, [a, b, num::relative_error(a,b)].iter().format(" "));
	}
	}
	let dt = system.dt(pressure, &state);
	state = system.step(/*relative_tolerance:*/ 1e-4, /*absolute_tolerance:*/ 1e-14, time_step, pressure, state);
	println!("{}", species.iter().map(|k| format!("{:^15}",k)).format(" "));
	let dtn : &[_;S-1] = iter::vec::Suffix::suffix(&dt);
	for &net_production_rates in &[&iter::vec::eval(dtn, |dtn| dtn/volume), other_net_production_rates] { println!("{}", net_production_rates.iter().map(|c| format!("{:15.6e}",c)).format(" ")); }
	println!("{}", dtn.iter().map(|dtn| dtn/volume).zip(other_net_production_rates).format_with(" ", |(a,&b), f| f(&format_args!("{:15.9e}", a-b))));
	println!("{}", dtn.iter().map(|dtn| dtn/volume).zip(other_net_production_rates).format_with(" ", |(a,&b), f| f(&format_args!("{:15.9e}", num::relative_error(a,b)))));
	let (_, _, ref concentrations) : (_,_,[_;S]) = System::state(pressure, &state);
	println!("{}", &[concentrations, other_concentrations].iter().format_with("\n",|concentrations, f| f(&format_args!("{}",concentrations.iter().format_with(" ", |c, f| f(&format_args!("{:15.6e}",c)))))));
	println!("{}", concentrations.iter().zip(other_concentrations).format_with(" ", |(c,o), f| f(&format_args!("{:15.8e}", c-o))));
	println!("{}", concentrations.iter().zip(other_concentrations).format_with(" ", |(&c,&o), f| f(&format_args!("{:15.9e}", num::relative_error(c,o)))));
	println!("{:e}", concentrations.iter().zip(other_concentrations).map(|(&c,&o)| num::abs(c-o)/o).max_by(|a,b| PartialOrd::partial_cmp(a,b).unwrap()).unwrap());
	}
	println!("OK");*/
}
