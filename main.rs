#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds,iter_partition_in_place)]#![allow(non_snake_case)]
mod yaml; mod device;
use {iter::map, ast::*};
fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(std::env::args().skip(1).next().unwrap())?)?)?;
	let model = yaml::parse(&model);
	use crate::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	let reactions = map(&*model.reactions, |r| Reaction::new(species_names, r));
	let ref state = initial_state(&model);
	use itertools::Itertools;
	if true {
		let rates = reaction::rates(&species.thermodynamics, &reactions);
		let rates = device::assemble(&rates);
		assert!(state.volume == 1.);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let states_len = 65536*8;
		let active_amounts = map(&amounts[0..amounts.len()-1], |&n| vec![n as f32; states_len].into_boxed_slice());
		let_!{ [energy_rate_RT, rates @ ..] = &*rates(&[*pressure_R as f32], &*([&[&vec![total_amount as f32; states_len] as &[_], &vec![*T as f32; states_len] as &[_]] as &[_], &*map(&*active_amounts, |a| &**a)].concat()))? => {
		let all_same = |array:&[f32]| { assert!(array.len() == states_len); for &v in array { assert_eq!(v, array[0]); } array[0] };
		let dtT_T = all_same(energy_rate_RT) as f64;
		let rates = rates.iter().map(|a| all_same(a));
		eprintln!("{}, dtT_T: {dtT_T:.3e}", rates.zip(&**species_names).format_with(", ", |(rate, name), f| f(&format!("{name}: {rate:.0}").to_string())));
		}}
	}
	#[cfg(feature="transport")] {
		let transport = wrap(device, transport::properties::<4>(&species));
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let_!{ [viscosity, thermal_conductivity, mixture_diffusion_coefficients @ ..] = &*transport([*pressure_R], &[&[total_amount, *T], active_amounts].concat()) => {
			eprintln!("μ: {viscosity:.4e}, λ: {thermal_conductivity:.4}, D: {:.4e}", mixture_diffusion_coefficients.iter().format(" "));
		}}
	}
	Ok(())
}
