#![feature(format_args_capture,array_map,in_band_lifetimes,default_free_fn,associated_type_bounds)]#![allow(non_snake_case,non_upper_case_globals)]
mod yaml; mod device;
use {anyhow::Result, itertools::Itertools, std::env::*, device::*};
fn main() -> Result<()> {
	let path = args().skip(1).next().unwrap();
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path)?)?)?;
	let model = yaml::parse(&model);
	use combustion::*;
	let (ref species_names, ref species, _, reactions, ref state) = new(&model);

	if false {
		let rates = reaction::rates(&species.thermodynamics, &reactions);
		let rates = with_repetitive_input(assemble(&rates), 1<<19);
		assert!(state.volume == 1.);
		let State{temperature, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let_!{ [dtT_T, rates @ ..] = &*rates(&[*pressure_R], &[&[total_amount, *temperature], &amounts[0..amounts.len()-1]].concat())? => {
		eprintln!("{}, dtT_T: {dtT_T:.3e}", rates.iter().zip(&**species_names).format_with(", ", |(rate, name), f| f(&format!("{name}: {rate:.0}").to_string())));
		}}
	}
	#[cfg(feature="transport")] {
		let transport = transport::properties::<4>(&species);
		let transport = with_repetitive_input(assemble(&transport), 1);
		let State{temperature, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let_!{ [viscosity, thermal_conductivity, mixture_diffusion_coefficients @ ..] = &*transport(&[*pressure_R], &[&[total_amount, *temperature], &amounts[0..amounts.len()-1]].concat())? => {
			eprintln!("μ: {viscosity:.4e}, λ: {thermal_conductivity:.4}, D: {:.4e}", mixture_diffusion_coefficients.iter().format(" "));
		}}
	}
	Ok(())
}
