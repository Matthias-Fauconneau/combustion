#![feature(box_syntax)]

#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	use combustion::*;
	let Simulation{species, system, mut state} = Simulation::new(&system)?;
	let mut app = ui::app::App::new(ui::plot::Plot::new(box [&["T"] as &[_], &species], vec!(state.clone().into())))?;
	app.idle = box move |plot| {
		for _ in 0..100000 { state.step(&system); }
		plot.values.push(state.clone().into());
		let density = system.average_molar_mass * system.pressure / (ideal_gas_constant * state.temperature);
		use iter::from_iter;
		#[allow(non_snake_case)] let B = from_iter(system.thermodynamics.iter().map(|thermodynamic| thermodynamic.b(state.temperature)));
		let concentrations = from_iter(scale(density, mul(recip(system.molar_masses.iter().copied()), state.mass_fractions.iter().copied())));
		for reaction in system.reactions.iter() {
			let Rate{forward_base_rate, reverse_base_rate, efficiency} = reaction.rate(state.temperature, &B, &concentrations);
			let base_rate = forward_base_rate - reverse_base_rate;
			let net_rate = efficiency * base_rate;
			//println!("{:e}", net_rate);
			if net_rate > 1e-20 { /*println!("{:e}", net_rate);*/ return true; }
		}
		false
	};
	app.run()?
}
