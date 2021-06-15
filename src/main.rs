#![feature(array_methods, array_map, associated_type_bounds)]#![allow(non_snake_case, uncommon_codepoints)]
use {std::borrow::Borrow, itertools::Itertools};
pub fn pretty(iter: impl IntoIterator<Item:Borrow<f64>>) -> String {
	iter.into_iter().format_with(", ", |e, f| f(&format!("{:6.6}", float_pretty_print::PrettyPrintFloat(*e.borrow())))).to_string()
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = include_bytes!("../LiDryer.ron");
	use combustion::*;
	let model = model::Model::new(model)?;
	let (ref species_names, ref species) = Species::new(&model.species);
	let state = initial_state(&model);
	#[cfg(feature="transport")] {
		let ref transport_polynomials = species.transport_polynomials();
		let transport = transport::transport(&species.molar_mass, transport_polynomials, &state);
		println!("{:.1e}", transport.viscosity);
		println!("{:.2}", transport.thermal_conductivity);
		println!("{:.2}", transport.mixture_diffusion_coefficients.iter().format(", "));
	}
	#[cfg(feature="reaction")] {
		use reaction::*;
		let exp_Gibbs_RT = exp_Gibbs_RT(&species.thermodynamics[0..species.len()-1]);
		let rates = rates(&iter::map(&*model.reactions, |r| Reaction::new(species_names, r)));
		assert!(state.volume == 1.);
		let State{temperature: T, amounts: concentrations, ..} = state;
		let log_T = f64::log2(T);
		let T2 = T*T;
		let T3 = T*T2;
		let T4 = T*T3;
		let rcp_T = 1./T;
		let exp_Gibbs0_RT = exp_Gibbs_RT(&[log_T,T,T2,T3,T4,rcp_T],&[]);
		let P0_RT = NASA7::reference_pressure / T;
		eprintln!("{}", pretty(&*rates(&[log_T,T,T2,T4,rcp_T,num::sq(rcp_T),P0_RT,1./P0_RT], &[&exp_Gibbs0_RT, &concentrations])));
	}
	(species_names, species);
	Ok(())
}
