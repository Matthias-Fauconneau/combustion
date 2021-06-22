#![feature(array_methods, array_map, associated_type_bounds, format_args_capture)]#![allow(non_snake_case)]
pub fn dot(a: &[f64], b: &[f64]) -> f64 { assert!(a.len()==b.len()); iter::dot(a.iter().copied().zip(b.iter().copied())) }

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = yaml_model::Loader::load_from_str(std::str::from_utf8(&std::fs::read(std::env::args().skip(1).next().unwrap())?)?)?;
	let model = yaml_model::parse(&model)?;
	use chemistry::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	let ref state = initial_state(&model);
	use {iter::map, itertools::Itertools, ast::wrap};
	if true {
		use reaction::*;
		let rates = wrap(rates(species.thermodynamics, &iter::map(&*model.reactions, |r| Reaction::new(species_names, r))));
		assert!(state.volume == 1.);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let rates = rates([T, pressure_R],[amounts]);
		let (energy_rate_RT, rates) = (rates[0], rates[1..]);
		eprintln!("{}, HRR: {:.3e}", rates.iter().zip(&**species_names).format_with(", ", |(rate, name), f| f(&format!("{name}: {rate:.0}").to_string())), NA * kB * T * -energy_rate_RT);
	}
	#[cfg(feature="transport")] {
		use transport::*;
		let ref transport = Polynomials::<4>::new(&species);
		let viscosity_T_12 = wrap(viscosity_T_12(&species.molar_mass, &transport.sqrt_viscosity_T_14));
		let thermal_conductivity_T_12_2 = wrap(thermal_conductivity_T_12_2(&transport.thermal_conductivity_T_12));
		let P_T_32_mixture_diffusion_coefficients = wrap(P_T_32_mixture_diffusion_coefficients(&map(&*transport.binary_thermal_diffusion_coefficients_T32, |row| &**row)));
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let T_12 = f64::sqrt(*T);
		let ln_T = f64::ln(*T);
		let ln_T_2 = ln_T*ln_T;
		let ln_T_3 = ln_T_2*ln_T;
		let total_amount = amounts.iter().sum::<f64>();
		let mole_fractions = map(&**amounts, |n| n / total_amount);
		let molar_mass = dot(&mole_fractions, &species.molar_mass);
		let viscosity = T_12*viscosity_T_12([ln_T, ln_T_2, ln_T_3], [&mole_fractions])[0];
		let thermal_conductivity = (T_12/2.)*thermal_conductivity_T_12_2([ln_T, ln_T_2, ln_T_3], [&mole_fractions])[0];
		let mass_fractions =	map(mole_fractions.iter().zip(&*species.molar_mass), |(x,m)| x * m / molar_mass);
		let T_32_P = T*T_12/(pressure_R*NA*kB);
		let P_T_32_mixture_diffusion_coefficients = P_T_32_mixture_diffusion_coefficients([ln_T, ln_T_2, ln_T_3], [&mole_fractions, &mass_fractions]);
		let mixture_diffusion_coefficients = map(&*P_T_32_mixture_diffusion_coefficients, |P_T_32_D| T_32_P*P_T_32_D);
		eprintln!("μ: {viscosity:.4e}, λ: {thermal_conductivity:.4}, D: {:.4e}", mixture_diffusion_coefficients.iter().format(" "));
	}
	Ok(())
}
