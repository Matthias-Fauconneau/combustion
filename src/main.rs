#![feature(bindings_after_at)]

#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	use combustion::*;
	let Simulation{species: specie_names, system, mut state} = Simulation::new(&system)?;
	while state.time < 4e-6 {
		for _ in 0..10000 { state.step(&system); }
	}
	dbg!(state.temperature);
	/*use iter::from_iter;
	let density = system.average_molar_mass * system.pressure / (ideal_gas_constant * state.temperature);
	dbg!(density);
	let specific_enthalpy : f64 = mul(state.mass_fractions.iter().copied(), system.thermodynamics.iter().zip(system.molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_enthalpy(state.temperature) / molar_mass)).sum();
	dbg!(specific_enthalpy);
	let specific_entropy : f64 = mul(state.mass_fractions.iter().copied(), system.thermodynamics.iter().zip(system.molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_entropy(state.temperature) / molar_mass)).sum();
	dbg!(specific_entropy);
	let specific_heat_capacity : f64 = mul(state.mass_fractions.iter().copied(), system.thermodynamics.iter().zip(system.molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_heat_capacity(state.temperature) / molar_mass)).sum();
	dbg!(specific_heat_capacity);
	#[allow(non_snake_case)] let B = from_iter(system.thermodynamics.iter().map(|thermodynamic| thermodynamic.b(state.temperature)));
	let concentrations = from_iter(scale(density, mul(recip(system.molar_masses.iter().copied()), state.mass_fractions.iter().copied())));*/
	for _reaction@Reaction{equation,..} in system.reactions.iter() {
		use itertools::Itertools;
		let [left, right] = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|(species, coefficients)| species.iter().zip(coefficients.iter()).map(|(&specie, coefficient)| (specie_names[specie], coefficient)).format(" ")));
		/*let Rate{forward_base_rate, reverse_base_rate, efficiency} = reaction.rate(state.temperature, &B, &concentrations);
		let forward_rate = efficiency * forward_base_rate;
		let reverse_rate = efficiency * reverse_base_rate;
		let net_rate = forward_rate - reverse_rate;
		print!("{:14.5e}{:14.5e}{:14.5e} {:32}", forward_rate, reverse_rate, net_rate);*/
		println!("{:?} = {:?}", left, right);
	}
}
