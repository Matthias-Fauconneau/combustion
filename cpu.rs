#![feature(format_args_capture,iter_is_partitioned)]#![allow(non_snake_case)]

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = yaml_model::Loader::load_from_str(std::str::from_utf8(&std::fs::read(std::env::args().skip(1).next().unwrap())?)?)?;
	let model = yaml_model::parse(&model);
	use chemistry::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	let reactions = map(&*model.reactions, |r| Reaction::new(species_names, r));
	let active = {
		let active = map(0..species.len()-1, |k| reactions.iter().any(|Reaction{net,..}| net[k] != 0));
		assert!(active.iter().is_partitioned(|&active| active),
								"Species must be partionned with all active species first and all inert species last so that the state only consists of active species when solving kinetics: {active:?}");
		active.iter().position(|active| !active).unwrap()
	};
	let ref state = initial_state(&model);
	use {iter::map, itertools::Itertools, ast::{wrap, let_}};
	/*if true {
		let reaction_rates = wrap(reaction::reaction_rates_function(&species.thermodynamics[..species.len()-1], &map(&*model.reactions, |r| Reaction::new(species_names, r))));
		assert!(state.volume == 1.);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		assert!(*T == 1000. && *pressure_R == 101325./(kB*NA),  "{T} {pressure_R}");
		let total_amount = amounts.iter().sum::<f64>();
		let active_amounts = &amounts[0..amounts.len()-1];
		let reaction_rates = &*reaction_rates(&[&[*pressure_R, total_amount, *T], active_amounts].concat());
		println!("{:.0}", reaction_rates.iter().format(" "))
	}*/
	if true {
		let rates = wrap(reaction::rates(species.len(), &species.thermodynamics[0..active], &reactions));
		assert!(state.volume == 1.);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let active_amounts = &amounts[0..amounts.len()-1];
		let_!{ [energy_rate_RT, rates @ ..] = &*rates(&[&[*pressure_R, total_amount, *T], active_amounts].concat()) => {
		println!("{}, HRR: {:.3e}", rates.iter().zip(&**species_names).format_with(", ", |(rate, name), f| f(&format!("{name}: {rate:.0}").to_string())), NA * kB * T * -energy_rate_RT);
		}}
	}
	#[cfg(feature="transport")] {
		let transport = wrap(transport::properties::<4>(&species));
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let active_amounts = &amounts[0..amounts.len()-1];
		let_!{[viscosity, thermal_conductivity, mixture_diffusion_coefficients @ ..] = &*transport(&[&[*pressure_R,total_amount, *T], active_amounts].concat()) => {
		eprintln!("μ: {viscosity:.4e}, λ: {thermal_conductivity:.4}, D: {:.4e}", mixture_diffusion_coefficients.iter().format(" "));
		}}
	}
	Ok(())
}
