#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	use combustion::*;
	let Simulation{species: specie_names, system, mut state} = Simulation::new(&system)?;
	while state.time < 4e-6 {
		for _ in 0..10000 { state.step(&system); }
	}
	dbg!(state.temperature);
	let density = system.average_molar_mass * system.pressure / (ideal_gas_constant * state.temperature);
	dbg!(density);
	use iter::from_iter;
	let specific_enthalpy : f64 = mul(state.mass_fractions.iter(), system.thermodynamics.iter().zip(system.molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_enthalpy(state.temperature) / molar_mass)).sum();
	dbg!(specific_enthalpy);
	let specific_entropy : f64 = mul(state.mass_fractions.iter(), system.thermodynamics.iter().zip(system.molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_entropy(state.temperature) / molar_mass)).sum();
	dbg!(specific_entropy);
	let specific_heat_capacity : f64 = mul(state.mass_fractions.iter(), system.thermodynamics.iter().zip(system.molar_masses.iter()).map(|(thermodynamic, molar_mass)| thermodynamic.specific_heat_capacity(state.temperature) / molar_mass)).sum();
	dbg!(specific_heat_capacity);
	let concentrations = from_iter(scale(density, mul(recip(&system.molar_masses), state.mass_fractions.iter())));
	for Reaction{equation,model,..} in system.reactions.iter() {
		use itertools::Itertools;
		let [left, right] = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|(species, coefficients)| species.iter().zip(coefficients.iter()).map(|(&specie, coefficient)| (specie_names[specie], coefficient)).format(" ")));
		let [forward, reverse] = iter::array::Iterator::collect::<[_;2]>(equation.iter().map(|(species, coefficients)| rate((&species, &coefficients), model, state.temperature, &concentrations)));
		let net = forward - reverse;
		println!("{:14.5e}{:14.5e}{:14.5e} {:32}", forward, reverse, net, format!("{:?} = {:?}", left, right));
	}
}
