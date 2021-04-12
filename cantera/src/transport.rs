use super::*;

pub fn check(model: model::Model, Simulation{state, ..}: &Simulation) {
	let (species_names, species) = Species::new(model.species);
	let transport_polynomials = species.transport_polynomials();

	let len = species.len();
	let file = std::ffi::CStr::from_bytes_with_nul(b"gri30.yaml\0").unwrap().as_ptr();
	let name = std::ffi::CStr::from_bytes_with_nul(b"gri30\0").unwrap().as_ptr();
	let phase = unsafe{thermo_newFromFile(file, name)};
	//{
		let cantera_species_names = (0..len).map(|k| {
			let mut specie = [0; 8];
			unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
			unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
		}).collect::<Box<_>>();
		/*assert_eq!(cantera_species_name.iter().map(String::as_str).collect::<Box<_>>(), species_names);
	}*/
	assert!(unsafe{thermo_nSpecies(phase)} == len);
	assert!(state.amounts.len() == len && !state.amounts.iter().any(|&n| n<0.));
	let cantera_order = |o: &[f64]| (0..o.len()).map(|i| o[species_names.iter().position(|&s| s==cantera_species_names[i]).unwrap()]).collect::<Box<_>>();
	unsafe{thermo_setMoleFractions(phase, state.amounts.len(), cantera_order(&state.amounts).as_ptr(), 1)}; // /!\ Needs to be set before pressure
	dbg!(state.pressure * NA, state.temperature);
	unsafe{thermo_setTemperature(phase, state.temperature)};
	unsafe{thermo_setPressure(phase, state.pressure * NA)}; // /!\ Needs to be set after mole fractions
	let transport = unsafe{trans_newDefault(phase, 5)};
	let viscosity = unsafe{trans_viscosity(transport)};
	let thermal_conductivity  = unsafe{trans_thermalConductivity(transport)};
	//let ref binary_thermal_diffusion_coefficients = {let mut array = vec![0.; len*len]; unsafe{trans_getBinDiffCoeffs(transport, len as i32, array.as_mut_ptr())}; array};
	let ref mixture_mass_averaged_thermal_diffusion_coefficients = {
		let mut array = vec![0.; len];
		let iok = unsafe{trans_getThermalDiffCoeffs(transport, len as i32, array.as_mut_ptr())};
		assert!(iok == 0, "{}", iok);
		if false {
			array.into_boxed_slice()
		} else {
			if false { // N2
				"5.88830303e-04 9.74876528e-04 2.54222365e-04 1.48570420e-04 2.49489253e-04 2.11459327e-04 1.65222080e-04 1.64120812e-04 2.38722194e-04 2.71708310e-04 1.85500267e-04 1.85500267e-04 1.81398673e-04 1.97357482e-04 1.63943507e-04 1.32876379e-04 1.41676796e-04 1.40531852e-04 1.38019881e-04 1.38019881e-04 1.37553058e-04 1.39820441e-04 1.38417244e-04 1.37105272e-04 1.36668970e-04 1.25653611e-04 1.24639892e-04 2.05464735e-04 1.19597783e-04 1.19597783e-04 2.26754242e-04 2.67747169e-04 2.62329021e-04 1.99451789e-04 1.60259336e-04 1.62640124e-04 1.44415354e-04 1.31154960e-04 1.64801742e-04 1.61420566e-04 1.40439564e-04 1.39179409e-04 2.05460769e-04 1.31725863e-04 1.31725863e-04 1.31725863e-04 1.32333100e-04 1.64317769e-04 1.61220023e-04 9.72742481e-05 9.68451144e-05 1.19049476e-04 1.18523749e-04"
			} else {
				"6.01596776e-04 1.00226482e-03 2.55056423e-04 1.64626535e-04 2.50065057e-04 2.11551373e-04 1.62992833e-04 1.61816274e-04 2.39813147e-04 2.73446037e-04 1.85330651e-04 1.85330651e-04 1.81041369e-04 2.00619281e-04 1.62086393e-04 1.30003634e-04 1.39490186e-04 1.38275792e-04 1.35835729e-04 1.35835729e-04 1.35175678e-04 1.38074477e-04 1.36591137e-04 1.35202900e-04 1.34660395e-04 1.23590104e-04 1.22515786e-04 2.02913487e-04 1.16797862e-04 1.16797862e-04 2.27286765e-04 2.69056404e-04 2.63346607e-04 1.99152991e-04 1.58283213e-04 1.60618127e-04 1.41421583e-04 1.28295471e-04 1.62720202e-04 1.59699277e-04 1.38351491e-04 1.37018153e-04 2.02909207e-04 1.28907824e-04 1.28907824e-04 1.28907824e-04 1.29558757e-04 1.63492991e-04 1.36504022e-04 9.47262305e-05 9.42679458e-05 1.16211421e-04 1.15648787e-04"
			}
			.split(' ').map(|s| s.parse().unwrap()).collect::<Box<[f64]>>()
		}
	};
	let order = |o: &[f64]| (0..o.len()).map(|i| o[cantera_species_names.iter().position(|s| s==species_names[i]).unwrap()]).collect::<Box<_>>();
	let mixture_mass_averaged_thermal_diffusion_coefficients = order(&mixture_mass_averaged_thermal_diffusion_coefficients);

	let ref transport = combustion::transport::transport(&species.molar_mass, &transport_polynomials, state);
	dbg!(&transport);

	dbg!(viscosity);
	dbg!(transport.viscosity);
	use {num::{sq, sqrt}, iter::{dot, zip, into::IntoCopied}};
	let State{temperature, pressure, volume, amounts, ..} = state;
	let T = *temperature;
	let amount = pressure * volume / (K * T);
	assert!(dbg!(num::abs(amounts.iter().sum::<f64>()-amount)) < 2e-6);
	let mole_fractions = amounts.iter().map(|n| n/amount).collect::<Box<_>>();
	let direct_viscosity = dot(zip(mole_fractions.copied(), |k|
		species.viscosity(k, T) /
		dot(zip(mole_fractions.copied(), |j|
			sq(1. + sqrt(species.viscosity(k, T))/sqrt(species.viscosity(j, T)) * sqrt(sqrt(species.molar_mass[j]/species.molar_mass[k]))) /
			(sqrt(8.) * sqrt(1. + species.molar_mass[k]/species.molar_mass[j]))
		))
	));
	dbg!(direct_viscosity);
	dbg!(num::relative_error(transport.viscosity, viscosity));
	assert!(num::relative_error(transport.viscosity, viscosity) < 5e-7);

	dbg!(thermal_conductivity);
	dbg!(transport.thermal_conductivity);
	let direct_thermal_conductivity = 1./2. * (
		dot(zip(mole_fractions.copied(), |k| species.thermal_conductivity(k, T))) +
		1. / dot(zip(mole_fractions.copied(), |k| 1. / species.thermal_conductivity(k, T)))
	);
	dbg!(direct_thermal_conductivity);
	dbg!(num::relative_error(transport.thermal_conductivity, thermal_conductivity));
	assert!(num::relative_error(transport.thermal_conductivity, thermal_conductivity) < 0.024);

	/*for k in 0..len { for j in 0..len {
		let r = binary_thermal_diffusion_coefficients[len*j+k];
		let v = transport_polynomials.binary_thermal_diffusion_coefficient(k, j, state.temperature) / (pressure*NA);
		assert!(num::relative_error(r, v) < 0.42, "{}", num::relative_error(r, v));
	}}*/

	let mean_molecular_weight = species.molar_mass.iter().sum::<f64>() / species.molar_mass.len() as f64;
	let transport_mixture_mass_averaged_thermal_diffusion_coefficients = iter::eval(mole_fractions.len(), |k|
		(1. - mole_fractions[k]*species.molar_mass[k]/mean_molecular_weight) /
		dot(zip(mole_fractions.copied(), |j| if j != k { 1. / (transport_polynomials.binary_thermal_diffusion_coefficient(k, j, T) / (pressure*NA))  } else { 0. }))
	);

	//dbg!(mixture_mass_averaged_thermal_diffusion_coefficients);
	//dbg!(&transport_mixture_mass_averaged_thermal_diffusion_coefficients);
	for (i, (&a, &b)) in
		species_names.iter().zip(mixture_mass_averaged_thermal_diffusion_coefficients.iter().zip(transport_mixture_mass_averaged_thermal_diffusion_coefficients.iter())) {
		//dbg!(a, b, num::relative_error(a,b));
		//println!("{} {} {:e}", a, b, num::abs(a-b));
		assert!(num::relative_error(a, b) < 0.14, "{}: {}", i, num::relative_error(a, b));
		assert!(num::abs(a-b) < 6e-5, "{}: {:e}", i, num::abs(a-b));
	}
}
