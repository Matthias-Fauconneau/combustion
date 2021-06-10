#![feature(array_methods, array_map, associated_type_bounds)]#![allow(non_snake_case, uncommon_codepoints)]
use {std::borrow::Borrow, itertools::Itertools};
pub fn pretty(iter: impl IntoIterator<Item:Borrow<f64>>) -> String { iter.into_iter().format_with(", ", |e, f| f(&float_pretty_print::PrettyPrintFloat(*e.borrow()))).to_string() }
fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = include_bytes!("../LiDryer.ron");
	use combustion::*;
	let model = model::Model::new(model)?;
	let (ref species_names, ref species) = Species::new(&model.species);
	#[cfg(feature="transport")] {
		let ref transport_polynomials = species.transport_polynomials();
		println!("{}", transport_polynomials.thermal_conductivity_T12.iter().format_with(",\n", |p, f| f(&format_args!("[{}]", pretty(p)))));
		/*println!("{}", transport_polynomials.binary_thermal_diffusion_coefficients_T32.iter().format_with(",\n",
				|c, f| f(&format_args!("[\n{}\n]", c.iter().format_with(",\n", |p, f| f(&format_args!("[{}]", pretty(p)))))) ));*/
		//println!("{}", transport::transport(&species.molar_mass, transport_polynomials, state).mixture_diffusion_coefficients);*/
	}
	#[cfg(feature="reaction")] {
		use reaction::*;
		let exp_Gibbs_RT = exp_Gibbs_RT(&species.thermodynamics[0..species.len()-1]);
		let rates = rates(&iter::map(&*model.reactions, |r| Reaction::new(species_names, r)));
		let state = initial_state(&model);
		assert!(state.volume == 1.);
		let State{temperature: T, amounts: concentrations, ..} = state;
		let log_T = f64::log2(T);
		let T2 = T*T;
		let T3 = T*T2;
		let T4 = T*T3;
		let rcp_T = 1./T;
		let exp_Gibbs0_RT = exp_Gibbs_RT(&[log_T,T,T2,T3,T4,rcp_T],&[]);
		let P0_RT = NASA7::reference_pressure / T;
		println!("{}", pretty(rates(&[log_T,T,T2,T4,rcp_T,num::sq(rcp_T),P0_RT,1./P0_RT], &[&exp_Gibbs0_RT, &concentrations])));
	}
	(species_names, species);
	Ok(())
}
