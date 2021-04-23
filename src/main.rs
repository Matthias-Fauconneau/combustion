#![feature(type_ascription)]#![allow(confusable_idents,non_snake_case,unused_variables,unused_mut,uncommon_codepoints)]
use {fehler::throws, error::Error, combustion::*};

#[throws] fn main() {
	let model = &std::fs::read("CH4+O2.ron")?;
	let ref model = model::Model::new(&model)?;
	let ref state = combustion::initial_state(model);
	let (ref species_names, ref species) = Species::new(&model.species);
	#[cfg(feature="transport")] {
		let ref transport_polynomials = species.transport_polynomials();
		/*use itertools::Itertools;
		use transport::*;
		/*for ((((T⃰, Ω22), A), B), C) in header_T⃰.iter().zip(Ω⃰22.iter()).zip(A⃰.iter()).zip(B⃰.iter()).zip(C⃰.iter()).skip(1) {
			println!("delta* fit at T* = {}", T⃰);
			println!("omega22 = [{:.5}]", Ω22.iter().format(", "));
			println!("A* = [{:.5}]", A.iter().format(", "));
			println!("B* = [{:.5}]", B.iter().format(", "));
			println!("C* = [{:.5}]", C.iter().format(", "));
		}*/
		/*let TransportPolynomials{sqrt_viscosity_T14, ..} = &transport_polynomials;
		for (specie, sqrt_viscosity_T14) in species_names.iter().zip(sqrt_viscosity_T14.iter()) { println!("{} = [{:.5}]", specie, sqrt_viscosity_T14.iter().format(", ")); }*/*/
		let transport = |single_specie: &str, temperature_C| {
				let pressure = 1e5/(K*NA);
				let temperature = 273.15+temperature_C;
				let volume = 1.;
				let amount = pressure * volume / (temperature * K);
				transport::transport(&species.molar_mass, transport_polynomials,
												&State{volume, temperature, pressure,
												amounts: {assert!(species_names.contains(&single_specie)); species_names.iter().map(|specie| if specie==&single_specie {amount} else {0.}).collect()}
												//amounts: {let amounts=vec![0.; species_names.len()].into_boxed_slice(); amounts[species_names.iter().position(|&specie| specie==single_specie).unwrap()] = amount; amounts}
												} )
		};
		let viscosity = |single_specie: &str, T, expected| {
			let e = f64::abs(transport(single_specie, T).viscosity*1e6-expected)/expected;
			//println!("{} {} {} {}", single_specie, T, transport(single_specie, T).viscosity, e);
			assert!(e < 0.03);
		};
		viscosity("AR", 25., 22.58); // 23
		viscosity("CH4", 25., 11.07); // 11.4
		let thermal_conductivity = |single_specie, T, expected| {
			let value = transport(single_specie, T).thermal_conductivity;
			let e = f64::abs(value-expected)/expected;
			//println!("{}", e);
			assert!(e < 0.001, "{} {}", value, expected);
		};
		thermal_conductivity("AR",500., 0.03650);
	}
	#[cfg(feature="reaction")] {
		use reaction::*;
		let reactions = iter::map(&*model.reactions, |r| Reaction::new(species_names, r));
		println!("{}", rate::<_,{Property::Volume}>(species, &*reactions)?);
		/*let (_, rate) = rate(species, &*reactions);
		let mut derivative = /*Derivative*/StateVector::<{Property::Volume}>(std::iter::repeat(0.).take(2+species.len()-1).collect());

		{
			let mut debug = vec![f64::NAN; model.reactions.len()*2].into_boxed_slice();
			let rate = {rate(state.constant(), &state.into(), &mut derivative, &mut debug); &derivative.0};
			let (_dt_temperature, _dt_variable, dt_amounts) =
				if let [_dt_temperature, _dt_variable, dt_amounts @..] = &rate[..] { (_dt_temperature, _dt_variable, dt_amounts) } else { unreachable!() };
			let ref molar_rate = dt_amounts.iter().map(|dtn| dtn/state.volume).collect::<Box<_>>();
			for (i, (specie, rate)) in species_names.iter().zip(molar_rate.iter()).enumerate() { if rate.abs() > 1e-29 { println!("{:3} {:4} {:>+0.3e}", i, specie, rate); } }
			println!("");
		}
		/*let len = 100000;
		let constant = state.constant();
		let state = state.into();
		let start = std::time::Instant::now();
		for _ in 0..len { rate(constant, &state, &mut derivative) }
		let end = std::time::Instant::now();
		let time = (end-start).as_secs_f32();
		println!("{:.1}ms\t{:.0}K/s", time*1e3, (len as f32)/time/1e3);*/*/
	}
}
