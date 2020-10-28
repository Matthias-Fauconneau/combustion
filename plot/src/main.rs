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
		let concentrations = from_iter(scale(density, mul(recip(system.molar_masses.iter().copied()), state.mass_fractions.iter().copied())));
		for Reaction{equation,model,..} in system.reactions.iter() {
			let [forward, reverse] = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|(species, coefficients)| rate((&species, &coefficients), model, state.temperature, &concentrations)));
			let net = forward - reverse;
			if net > 1e-5/*13*/ { println!("{:e}", net); return true; }
		}
		false
	};
	app.run()?
}
