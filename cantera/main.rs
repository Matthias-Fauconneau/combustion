#![allow(mixed_script_confusables, non_snake_case, incomplete_features, confusable_idents)]
#![feature(type_ascription, array_map, non_ascii_idents, const_generics, const_evaluatable_checked)]
use {fehler::throws, anyhow::Error};
use combustion::*;

#[throws] fn main() {
	let system = std::fs::read("CH4+O2.ron")?;
	const S: usize  = 53;
	test_reaction_cantera(&Simulation::<S>::new(&system)?)?;
	//test_transport_cantera(&Simulation::<S>::new(&system)?)?;
}

mod cantera {
extern "C" {
pub fn reaction(pressure: &mut f64, temperature: &mut f64, mole_proportions: *const std::os::raw::c_char, time_step: f64,
													species_len: &mut usize,
														species: &mut *const *const std::os::raw::c_char,
														net_productions_rates: &mut *const f64,
													reactions_len: &mut usize,
														equations: &mut *const *const std::os::raw::c_char,
														equilibrium_constants: &mut *const f64,
														forward: &mut *const f64,
														reverse: &mut *const f64,
													concentrations: &mut *const f64); }
extern "C" {
pub fn transport(pressure: f64, temperature: f64, mole_proportions: *const std::os::raw::c_char,
														viscosity: &mut f64, thermal_conductivity: &mut f64,
														species_len: &mut usize, species: &mut *const *const std::os::raw::c_char, mixture_averaged_thermal_diffusion_coefficients: &mut *const f64); }
}

#[throws] fn test_reaction_cantera<const S: usize>(Simulation{species_names, system, pressure_R, time_step, state, ..}: &Simulation<S>) where [(); S-1]:, [(); 1+S-1]: {
	assert!(System::<S>::volume == 1.);
	use {iter::Suffix, std::convert::TryInto, itertools::Itertools, num::relative_error};
	let (equations, equilibrium_constants, [forward, reverse], ref cantera_rate, ref _cantera_concentrations) = {
		let mole_proportions = format!("{}", species_names.iter().zip(&state.amounts).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
		let mole_proportions = std::ffi::CString::new(mole_proportions)?;
		use std::ptr::null;
		let mut pressure = pressure_R * (combustion::kB*combustion::NA);
		let mut temperature = state.temperature;
		//let (mut species_len, mut specie_names, mut net_productions_rates, mut concentrations) = (0, null(), null(), null());
		let (mut species_len, mut specie_names, mut net_productions_rates, mut concentrations, mut reactions_len, mut equations, mut equilibrium_constants, [mut forward, mut reverse])
			= (0, null(), null(), null(), 0, null(), null(), [null(); 2]);
		unsafe {
			cantera::reaction(&mut pressure, &mut temperature, mole_proportions.as_ptr(), *time_step, &mut species_len, &mut specie_names, &mut net_productions_rates,
																	&mut reactions_len, &mut equations, &mut equilibrium_constants, &mut forward, &mut reverse,
																	&mut concentrations);
			let specie_names = iter::box_collect(std::slice::from_raw_parts(specie_names, species_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
			let order = |o:&[_]| iter::vec::eval(species_names, |s| o[specie_names.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}", s, species_names))]);
			let net_productions_rates = order(std::slice::from_raw_parts(net_productions_rates, species_len)).map(|c| c*1000.); // kmol/m^3/s => mol/s [1m^3]
			let concentrations = order(std::slice::from_raw_parts(concentrations, species_len)).map(|c| c*1000.); // kmol/m^3 => mol/m^3
			(
				iter::box_collect(std::slice::from_raw_parts(equations, reactions_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap())),
				iter::box_collect(std::slice::from_raw_parts(equilibrium_constants, reactions_len).iter().zip(system.reactions.iter()).map(|(c,Reaction{Σnet, ..})| c*f64::powf(1e3, *Σnet))),
				[forward, reverse].map(|r| iter::box_collect(std::slice::from_raw_parts(r, reactions_len).iter().map(|c| c*1000.))),
				iter::vec::eval(net_productions_rates, |dtC| dtC * System::<S>::volume)[0..S-1].try_into()?:[_;S-1],
				//State{temperature, amounts: iter::vec::eval(concentrations, |c| c * System::<S>::volume)}
				concentrations
			)
		}
	};
	if true {
    let other_reactions = iter::box_collect(iter::zip!(equations, equilibrium_constants, forward, reverse).into_iter());
    let reactions = {
        let System{species: Species{thermodynamics, ..}, reactions, ..} = &system;
        let State{temperature, amounts} = state;
        use num::log;
        let T = *temperature;
        let logP0_RT = log(NASA7::reference_pressure_R) - log(T);
        use iter::{Prefix, vec::eval, eval};
        let ref H = eval(thermodynamics, |s| s.specific_enthalpy(T));
				let ref H_T = eval(H.prefix(), |H| H/T);
        let ref G = eval!(thermodynamics.prefix(), H_T; |s, h_T| h_T - s.specific_entropy(T)); // (H-TS)/RT
        let ref concentrations = eval(amounts, |&n| n / System::<S>::volume);
        let ref log_concentrations = eval(concentrations, |&c| f64::ln(c));
        iter::box_collect(reactions.iter().map(move |Reaction{reactants, products, rate_constant, model, net, Σnet, ..}| {
            let equation = format!("{}", [reactants, products].iter().format_with(" <=> ", |side, f| {
							//let other_species = ["Ar","H","H2","H2O","HO2","H2O2","O","O2","OH"];
							f(&side.iter().enumerate().filter(|(_,&ν)| ν > 0.)
								.sorted_by_key(|&(k,_)| species_names[k])
								.format_with(" + ", |(specie, &ν), f| if ν > 1. { f(&format_args!("{} {}",ν,species_names[specie])) } else { f(&species_names[specie]) }))?;
							use Model::*; match model {
									Elementary => {},
									ThreeBody{..} => f(&" + M")?,
									PressureModification{..}|Falloff{..} => f(&" (+M)")?,
							};
							Ok(())
            }));
            let log_kf = combustion::reaction::log_arrhenius(rate_constant, T);
						let c = model.efficiency(T, concentrations, log_kf);
            let mask = |mask, v| iter::zip!(mask, v).map(|(&mask, v):(_,&_)| if mask != 0. { *v } else { 0. });
						use iter::{into::IntoMap, vec::Dot};
						let Rf = f64::exp(reactants.dot(mask(reactants, log_concentrations)) + log_kf);
						let log_equilibrium_constant = -net.dot(G) + Σnet*logP0_RT;
						let Rr = f64::exp(products.dot(mask(products, log_concentrations)) + log_kf - log_equilibrium_constant);
						(equation, f64::exp(log_equilibrium_constant), c*Rf, c*Rr)
        }))
    };
    for (a, b) in reactions.iter().zip(other_reactions.iter()) {
        let ref e = a.0;
        //let b = other_reactions.into_iter().find(|(k,_,_,_)| k==&e).expect(&format!("{} {}", &e, other_reactions.iter().map(|b| b.0).format(" ")));
        assert_eq!(e, &b.0.replace(" =>"," <=>"));
        assert!(num::relative_error(a.1, b.1) < 0.002, "{}", num::relative_error(a.1, b.1));
				let a = a.2-a.3;
        let b = b.2-b.3;
        if num::relative_error(a,b) != 0. { println!("{:30} {:15.2e}", e, [a, b, num::relative_error(a,b)].iter().format(" ")); }
    }
	}
	let (rate, /*jacobian*/) = system.rate/*and_jacobian*/(*pressure_R, &(*state).into()).unwrap();
	let rate: &[_; S-1] = rate.suffix();
	let table = species_names.iter().zip(rate.iter().zip(cantera_rate)).filter(|(_,(&a,&b))| a != 0. || b != 0.).map(|(&header,(&a,&b))| {
		fn to_string(v: f64) -> String { if v == 0. { "0".to_owned() } else { format!("{:.0e}", v) } }
		let column = [header.to_owned(), to_string(a), to_string(b), to_string(relative_error(a,b))];
		let width = column.iter().map(|s| s.len()).max().unwrap();
		(column, width)
	}).collect::<Box<_>>();
	fn print<const R: usize>(table: &[([String; R], usize)]) { for row in 0..R { println!("{}", table.iter().format_with(" ", |(c,width), f| f(&format_args!("{:width$}", c[row], width=width)))); } }
	print(&table);
	assert!(rate == cantera_rate);

	/*let ref u: [f64; 1+S-1] = (*state).into();
	let mut cvode = cvode::CVODE::new(move |u| system.rate/*and_jacobian*/(*pressure_R, u).map(|(rate, /*jacobian*/)| rate), u);
	let (time, u) = cvode.step(*time_step, &(*state).into());
	assert_eq!(time, *time_step);
	let state = State::<S>::new(state.amounts.iter().sum(), &u);
	dbg!(state, cantera_state);
	assert!(state == cantera_state);*/
}

#[allow(dead_code)] #[throws] fn test_transport_cantera<const S: usize>(Simulation{system, state, pressure_R, species_names, ..}: &Simulation<S>) where [(); S-1]: {
	let transport = system.transport(*pressure_R, &state);
	let cantera = {
		use itertools::Itertools;
		let mole_proportions = format!("{}", species_names.iter().zip(&state.amounts).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
		let mole_proportions = std::ffi::CString::new(mole_proportions)?;
		use std::ptr::null;
		let ([mut viscosity, mut thermal_conductivity], mut species_len, mut specie_names, mut mixture_averaged_thermal_diffusion_coefficients) = ([0.; 2], 0, null(), null());
		unsafe {
			let pressure = pressure_R * (combustion::kB*combustion::NA);
			cantera::transport(pressure, state.temperature, mole_proportions.as_ptr(), &mut viscosity, &mut thermal_conductivity, &mut species_len, &mut specie_names,
																		&mut mixture_averaged_thermal_diffusion_coefficients);
			let specie_names = iter::box_collect(std::slice::from_raw_parts(specie_names, species_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
			let order = |o:&[_]| iter::vec::eval(species_names, |s| o[specie_names.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}", s, species_names))]);
			let mixture_averaged_thermal_diffusion_coefficients = order(std::slice::from_raw_parts(mixture_averaged_thermal_diffusion_coefficients, species_len));
			Transport{viscosity, thermal_conductivity, mixture_averaged_thermal_diffusion_coefficients}
		}
	};
	dbg!(&transport, &cantera);
	dbg!((transport.viscosity-cantera.viscosity)/cantera.viscosity, (transport.thermal_conductivity-cantera.thermal_conductivity)/cantera.thermal_conductivity);
	assert!(f64::abs(transport.viscosity-cantera.viscosity)/cantera.viscosity < 0.03, "{}", transport.viscosity/cantera.viscosity);
	assert!(f64::abs(transport.thermal_conductivity-cantera.thermal_conductivity)/cantera.thermal_conductivity < 0.05, "{:?}", (transport.thermal_conductivity, cantera.thermal_conductivity));
}
