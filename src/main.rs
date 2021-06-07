use combustion::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = &std::fs::read("LiDryer.ron")?;
	let model = model::Model::new(&model)?;
	let ref state = combustion::initial_state(&model);
	let (ref species_names, ref species) = combustion::Species::new(&model.species);
	#[cfg(feature="transport")] {
		println!("{:#?}", &species);
		let ref transport_polynomials = species.transport_polynomials();
		//println!("{:#?}", &transport_polynomials.binary_thermal_diffusion_coefficients_T32);
		//println!("{}", transport::transport(&species.molar_mass, transport_polynomials, state).mixture_diffusion_coefficients);
	}
	#[cfg(all(feature="reaction", feature="jit"))] {
		use reaction::{*, Property::*};
		let reactions = iter::map(&*model.reactions, |r| Reaction::new(species_names, r));
		//rate::<_,_,{Property::Volume}>(species, &*reactions, &mut std::io::stdout())?
		let rate = rate_function(species, &*reactions);
		let mut derivative = /*Derivative*/StateVector::<{Volume}>(std::iter::repeat(0.).take(2+species.len()-1).collect());

		{
			let rate = {rate(state.constant(), &state.into(), &mut derivative); &derivative.0};
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
		println!("{:.1}ms\t{:.0}K/s", time*1e3, (len as f32)/time/1e3);*/
	}
	Ok(())
}
