#![feature(box_syntax)]

#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	use combustion::*;
	let simulation @ Simulation{species_names, time_step, state, ..} = Simulation::<35>::new(&system)?;
	let system = CVODE::new(&simulation);
	fn from<const S: usize>(State{temperature, amounts}: State<S>) -> Box<[Box<[f64]>]> { box [box [temperature] as Box<[_]>, amounts] as Box<[_]> }
	/*impl<const S: usize> From<State<S>> for [f64; 1+S-1] /*where [(); S-1]:, [(); 1+S-1]:*/ {
	fn from(State{temperature, amounts}: State) -> Self { use iter::{array_from_iter as from_iter, into::IntoChain}; from_iter([*temperature].chain(amounts[..S-1].try_into().unwrap():[_;S-1])) }
}*/
	//fn from*
	let mut app = ui::app::App::new(ui::plot::Plot::new(box [&["T"] as &[_], &species_names], vec![(0., from(state))]))?;
	let mut time = 0.;
	app.idle = box move |plot| {
		/*for _ in 0..100000*/ { system.step(time_step, &state.into()); time += system.time_step; }
		plot.values.push( (time*1e9/*ns*/, from(state)) );
		/*let density = system.average_molar_mass * system.pressure / (ideal_gas_constant * state.temperature);
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
		false*/
		true
	};
	app.run()?
}
