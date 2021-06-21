#![feature(array_methods, array_map, associated_type_bounds, format_args_capture)]#![allow(non_snake_case)]
pub fn dot(a: &[f64], b: &[f64]) -> f64 { assert!(a.len()==b.len()); iter::dot(a.iter().copied().zip(b.iter().copied())) }
use combustion::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = include_bytes!("../LiDryer.ron");
	let model = model::Model::new(model)?;
	let (ref species_names, ref species) = Species::new(&model.species);
	let ref state = initial_state(&model);
	#[cfg(feature="itertools")] if true {
		use {iter::map, itertools::Itertools, program::wrap, reaction::*};
		let exp_Gibbs_RT = wrap(exp_Gibbs_RT(&species.thermodynamics[0..species.len()-1]));
		let enthalpy_RT = wrap(enthalpy_RT(&species.thermodynamics[0..species.len()-1]));
		let rates = wrap(rates(&iter::map(&*model.reactions, |r| Reaction::new(species_names, r))));
		assert!(state.volume == 1.);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let log_T = f64::log2(*T);
		let T2 = T*T;
		let T3 = T*T2;
		let T4 = T*T3;
		let rcp_T = 1./T;
		let exp_Gibbs0_RT = exp_Gibbs_RT([log_T,*T,T2,T3,T4,rcp_T],[]);
		let P0_RT = NASA7::reference_pressure / T;
		let total_amount = amounts.iter().sum::<f64>();
		let total_concentration = pressure_R / T;
		let concentrations = map(&**amounts, |n| n/total_amount*total_concentration);
		let rates = rates([log_T,*T,T2,T4,rcp_T,num::sq(rcp_T),P0_RT,1./P0_RT], [&exp_Gibbs0_RT, &concentrations]);
		let enthalpy_RT = enthalpy_RT([log_T,*T,T2,T3,T4,rcp_T],[]);
		let energy_rate_RT = dot(&rates, &enthalpy_RT);
		eprintln!("{}, HRR: {:.3e}", rates.iter().zip(&**species_names).format_with(", ", |(rate, name), f| f(&format!("{name}: {rate:.0}").to_string())), NA * kB * T * -energy_rate_RT);
	}
	#[cfg(all(feature="transport",feature="itertools"))] {
		use {iter::map, itertools::Itertools, program::wrap, transport::*};
		let ref transport_polynomials = species.transport_polynomials();
		let viscosity_T_12 = wrap(viscosity_T_12(&species.molar_mass, &transport_polynomials.sqrt_viscosity_T_14));
		let thermal_conductivity_T_12_2 = wrap(thermal_conductivity_T_12_2(&transport_polynomials.thermal_conductivity_T_12));
		let P_T_32_mixture_diffusion_coefficients = P_T_32_mixture_diffusion_coefficients(&map(&*transport_polynomials.binary_thermal_diffusion_coefficients_T32, |row| &**row));
		#[cfg(feature="debug")] if false { println!("{:?}", P_T_32_mixture_diffusion_coefficients); }
		let P_T_32_mixture_diffusion_coefficients = wrap(P_T_32_mixture_diffusion_coefficients);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let T_12 = f64::sqrt(*T);
		let ln_T = f64::ln(*T);
		let ln_T_2 = ln_T*ln_T;
		let ln_T_3 = ln_T_2*ln_T;
		let total_amount = amounts.iter().sum::<f64>();
		let mole_fractions = iter::map(&**amounts, |n| n / total_amount);
		let molar_mass = dot(&mole_fractions, &species.molar_mass);
		let viscosity = T_12*viscosity_T_12([ln_T, ln_T_2, ln_T_3], [&mole_fractions])[0];
		let thermal_conductivity = (T_12/2.)*thermal_conductivity_T_12_2([ln_T, ln_T_2, ln_T_3], [&mole_fractions])[0];
		let mass_fractions = iter::map(mole_fractions.iter().zip(species.molar_mass.iter()), |(x,m)| x * m / molar_mass);
		let T_32_P = T*T_12/(pressure_R*NA*kB);
		let P_T_32_mixture_diffusion_coefficients = P_T_32_mixture_diffusion_coefficients([ln_T, ln_T_2, ln_T_3], [&mole_fractions, &mass_fractions]);
		let mixture_diffusion_coefficients = map(P_T_32_mixture_diffusion_coefficients, |P_T_32_D| T_32_P*P_T_32_D);
		eprintln!("μ: {viscosity:.4e}, λ: {thermal_conductivity:.4}, D: {:.4e}", mixture_diffusion_coefficients.iter().format(" "));
	}
	Ok(())
}
