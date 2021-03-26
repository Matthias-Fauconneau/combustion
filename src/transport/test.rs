use {super::{*, super::*}};

#[test] #[fehler::throws(error::Error)] fn test() {
	let model = &std::fs::read("CH4+O2.ron")?;
	let model = model::Model::new(&model)?;
	let ref state = Simulation::new(&model)?.state;
	let (species_names, species) = Species::new(model.species);
	let transport_polynomials = species.transport_polynomials();
	let transport = |single_specie, temperature_C| {
			let pressure = 1e5/(K*NA);
			let temperature = 273.15+temperature_C;
			let volume = 1.;
			let amount = pressure / temperature * volume;
			transport(&species.molar_mass, &transport_polynomials,
											&State{volume, temperature, pressure, amounts: species_names.iter().map(|&specie| if specie==single_specie {amount} else {0.}).collect()} )
	};
	let viscosity = |single_specie, T, expected| { let e = f64::abs(transport(single_specie, T).viscosity*1e6-expected)/expected; println!("{}", e); assert!(e < 0.07); };
	viscosity("Ar",25., 22.58);
	viscosity("CH4",25., 11.07);
	let thermal_conductivity = |single_specie, T, expected| { let got = transport(single_specie, T).thermal_conductivity; let e = f64::abs(got-expected)/expected; println!("{}", e);  assert!(e < 0.05, "{} {}", got, expected); };
	thermal_conductivity("Ar",500., 0.03650);
}
