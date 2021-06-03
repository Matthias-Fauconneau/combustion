#![feature(type_ascription)]#![allow(confusable_idents,non_snake_case,unused_variables,unused_mut,uncommon_codepoints)]
use combustion::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = &std::fs::read("H2.ron")?;
	let ref model = model::Model::new(&model)?;
	let ref state = combustion::initial_state(model);
	let (ref species_names, ref species) = Species::new(&model.species);
	#[cfg(feature="transport")] {
		let ref transport_polynomials = species.transport_polynomials();
		let pressure_R = 1e5/(K*NA);
		let temperature = 1000.;
		let volume = 1.;
		let amount = pressure_R * volume / temperature;
		println!("{}",
			transport::transport(&species.molar_mass, transport_polynomials, &State{volume, temperature, pressure_R, amounts: vec![1./(species.len() as f64); species.len()]})
			.mixture_diffusion_coefficients);
	}
	#[cfg(feature="reaction")] {
		use reaction::*;
		let reactions = iter::map(&*model.reactions, |r| Reaction::new(species_names, r));
		rate::<_,_,{Property::Volume}>(species, &*reactions, &mut std::io::stdout())?
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
	Ok(())
}
