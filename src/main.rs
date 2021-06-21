#![feature(array_methods, array_map, associated_type_bounds)]#![allow(non_snake_case, uncommon_codepoints)]
#[cfg(feature="float-pretty-print")] pub fn pretty(iter: impl IntoIterator<Item:std::borrow::Borrow<f64>>) -> String {
	use {std::borrow::Borrow, itertools::Itertools};
	iter.into_iter().format_with(", ", |e, f| f(&format!("{:6.6}", float_pretty_print::PrettyPrintFloat(*e.borrow())))).to_string()
}

use combustion::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = include_bytes!("../LiDryer.ron");
	let model = model::Model::new(model)?;
	let (ref species_names, ref species) = Species::new(&model.species); 	(species_names, species);
	let ref state = initial_state(&model);
	#[cfg(feature="float-pretty-print")] if false {
		use {program::wrap, reaction::*};
		let exp_Gibbs_RT = wrap(exp_Gibbs_RT(&species.thermodynamics[0..species.len()-1]));
		let rates = wrap(rates(&iter::map(&*model.reactions, |r| Reaction::new(species_names, r))));
		assert!(state.volume == 1.);
		let State{temperature: T, amounts: concentrations, ..} = state;
		let log_T = f64::log2(*T);
		let T2 = T*T;
		let T3 = T*T2;
		let T4 = T*T3;
		let rcp_T = 1./T;
		let exp_Gibbs0_RT = exp_Gibbs_RT([log_T,*T,T2,T3,T4,rcp_T],[]);
		let P0_RT = NASA7::reference_pressure / T;
		println!("{}", pretty(&*rates([log_T,*T,T2,T4,rcp_T,num::sq(rcp_T),P0_RT,1./P0_RT], [&exp_Gibbs0_RT, concentrations])));
	}
	#[cfg(all(feature="transport",feature="itertools"))] {
		use {iter::map, itertools::Itertools, program::wrap, transport::*};
		let ref transport_polynomials = species.transport_polynomials();
		let viscosity_T_12 = wrap(viscosity_T_12(&species.molar_mass, &transport_polynomials.sqrt_viscosity_T_14));
		let thermal_conductivity_T_12_2 = wrap(thermal_conductivity_T_12_2(&transport_polynomials.thermal_conductivity_T_12));
		let P_T_32_mixture_diffusion_coefficients = P_T_32_mixture_diffusion_coefficients(&map(&*transport_polynomials.binary_thermal_diffusion_coefficients_T32, |row| &**row));
		#[cfg(feature="debug")] println!("{:?}", P_T_32_mixture_diffusion_coefficients);
		let P_T_32_mixture_diffusion_coefficients = wrap(P_T_32_mixture_diffusion_coefficients);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let T_12 = f64::sqrt(*T);
		let ln_T = f64::ln(*T);
		let ln_T_2 = ln_T*ln_T;
		let ln_T_3 = ln_T_2*ln_T;
		let total_amount = amounts.iter().sum::<f64>();
		let mole_fractions = iter::map(&**amounts, |n| n / total_amount);
		let dot = |a:&[_],b:&[_]| iter::dot(a.iter().copied().zip(b.iter().copied()));
		let molar_mass = dot(&mole_fractions, &species.molar_mass);
		let viscosity = T_12*viscosity_T_12([ln_T, ln_T_2, ln_T_3], [&mole_fractions])[0];
		let thermal_conductivity = (T_12/2.)*thermal_conductivity_T_12_2([ln_T, ln_T_2, ln_T_3], [&mole_fractions])[0];
		let mass_fractions = iter::map(mole_fractions.iter().zip(species.molar_mass.iter()), |(x,m)| x * m / molar_mass);
		let T_32_P = T*T_12/(pressure_R*NA*kB);
		let P_T_32_mixture_diffusion_coefficients = P_T_32_mixture_diffusion_coefficients([ln_T, ln_T_2, ln_T_3], [&mole_fractions, &mass_fractions]);
		let mixture_diffusion_coefficients = map(P_T_32_mixture_diffusion_coefficients, |P_T_32_D| T_32_P*P_T_32_D);
		eprintln!("μ: {:.4e}, λ: {:.4}, D: {:.4e}", viscosity, thermal_conductivity, mixture_diffusion_coefficients.iter().format(" "));
	}
	Ok(())
}
