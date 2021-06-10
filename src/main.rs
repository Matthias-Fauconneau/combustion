#![feature(array_methods, array_map)]#![allow(non_snake_case, uncommon_codepoints)]
use combustion::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
	//let model = &std::fs::read("LiDryer.ron")?;
	let model = include_bytes!("../LiDryer.ron");
	let model = model::Model::new(model)?;
	let (ref _species_names, ref species) = Species::new(&model.species);
	use itertools::Itertools;
	fn pretty(iter: impl IntoIterator<Item=f64>) -> String { iter.into_iter().format_with(", ", |e, f| f(&float_pretty_print::PrettyPrintFloat(e))).to_string() }
	//let dbg = |label: &str, array:&[f64]| eprintln!("{}: {}", label, pretty(array));
	/*dbg("molar_mass", &species.molar_mass);
	dbg("internal_degrees_of_freedom", &species.internal_degrees_of_freedom);
	dbg("heat_capacity_ratio", &species.heat_capacity_ratio);
	dbg("well_depth_J", &species.well_depth_J);
	dbg("diameter", &species.diameter);
	dbg("permanent_dipole_moment", &species.permanent_dipole_moment);
	dbg("polarizability", &species.polarizability);
	dbg("rotational_relaxation", &species.rotational_relaxation);
	eprintln!("thermodynamics: {}", species.thermodynamics.iter().format_with(", ",
		|&NASA7{temperature_split, pieces}, f| f(&format_args!("({}, [{}])", temperature_split, pieces.iter().format_with(", ", |p, f| f(&format_args!("[{}]", pretty(p))))))
	));*/

	#[cfg(feature="transport")] {
		/*const*/let [temperature_min, temperature_max] : [f64; 2] = [300., 3000.];
		const N : usize = /*D+2 FIXME: Remez*/50;
		use iter::{into::{IntoCopied, IntoMap}, ConstSizeIterator};
		let T : [_; N] = iter::ConstRange.map(|n| temperature_min + (n as f64)/((N-1) as f64)*(temperature_max - temperature_min)).collect();
		//eprintln!("{}", pretty(T.iter().map(|&T| species.T⃰(0,0, T))));
		let (a, b, T) = (0, 0, T[47]);
		let ln_T⃰ = num::ln(species.T⃰ (a, b, T));
		/*const*/let header_ln_T⃰ = transport::header_T⃰.each_ref().map(|&T| num::ln(T));
		let interpolation_start_index = std::cmp::min((1+header_ln_T⃰ [1..header_ln_T⃰.len()].iter().position(|&header_ln_T⃰ | ln_T⃰ < header_ln_T⃰ ).unwrap())-1, header_ln_T⃰.len()-3);
		//dbg!(interpolation_start_index, header_ln_T⃰.len());
		if false {
			let ref transport_polynomials = species.transport_polynomials();
			eprintln!("{}", transport_polynomials.binary_thermal_diffusion_coefficients_T32.iter().format_with(",\n",
				|c, f| f(&format_args!("[\n{}\n]", c.iter().format_with(",\n", |p, f| f(&format_args!("[{}]", pretty(p.copied())))))) ));
			//println!("{}", transport::transport(&species.molar_mass, transport_polynomials, state).mixture_diffusion_coefficients);
		}
	}

	#[cfg(feature="reaction")] {
		use reaction::*;
		//println!("fn molar_heat_capacity_at_constant_pressure_R{}", molar_heat_capacity_at_constant_pressure_R(&species.thermodynamics));
		//println!("fn enthalpy_RT{}", enthalpy_RT(&species.thermodynamics));
		let exp_Gibbs_RT = exp_Gibbs_RT(&species.thermodynamics[0..species.len()-1]);
		let rates = rates(&iter::map(&*model.reactions, |r| Reaction::new(_species_names, r)));
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
		let rates = rates(&[log_T,T,T2,T4,rcp_T,num::sq(rcp_T),P0_RT,1./P0_RT], &[&exp_Gibbs0_RT, &concentrations]);
		use itertools::Itertools;
		eprintln!("{:e}", rates.iter().format(", "));
		//for rate in rates[species.len()-1..].iter() { eprintln!("     {:>+0.3e}", rate); }
		//let rates = &rates[..species.len()-1];
		//for (specie, rate) in _species_names.iter().zip(rates.iter()) { eprintln!("{:4} {:>+0.3e}", specie, rate); }
	}

	#[cfg(all(feature="jit"))] {
		use reaction::{*, Property::*};
		let rate = rate_function(species, model.reactions.map(|r| Reaction::new(species_names, r)));
		let mut derivative = /*Derivative*/StateVector::<{Volume}>(std::iter::repeat(0.).take(2+species.len()-1).collect());
		let ref state = initial_state(&model);
		{
			let rate = {rate(state.constant(), &state.into(), &mut derivative); &derivative.0};
			let (_dt_temperature, _dt_variable, dt_amounts) =
				if let [_dt_temperature, _dt_variable, dt_amounts @..] = &rate[..] { (_dt_temperature, _dt_variable, dt_amounts) } else { unreachable!() };
			let ref rates = dt_amounts.iter().map(|dtn| dtn/state.volume).collect::<Box<_>>();
			for (i, (specie, rate)) in _species_names.iter().zip(rates.iter()).enumerate() { if rate.abs() > 1e-29 { println!("{:3} {:4} {:>+0.3e}", i, specie, rate); } }
			println!("");
		}
		/*let len = 100000;
		let constant = state.constant();
		let state = state.into();
		let start = std::time::Instant::now();
		for _ in 0..len { rate(constant, &state, &mut derivative) }
		let end = std::time::Instant::now();
		let time = (end-start).as_secs_f32();
		println!("{:.1}ms\t{:.0}K/s", time*1e3, (len as f32)/time/1e3);*/
	}
	Ok(())
}
