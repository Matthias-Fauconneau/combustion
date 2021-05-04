use std::os::raw::c_char;
#[link(name = "cantera")]
extern "C" {
fn thermo_newFromFile(file_name: *const c_char, phase_name: *const c_char) -> i32;
fn thermo_nSpecies(n: i32) -> usize;
fn thermo_setTemperature(n: i32, t: f64) -> i32;
fn thermo_setMoleFractions(n: i32, len: usize, x: *const f64, norm: i32) -> i32;
fn thermo_getSpeciesName(n: i32, m: usize, len: usize, buffer: *mut c_char) -> i32;
fn thermo_setPressure(n: i32, p: f64) -> i32;
fn trans_newDefault(th: i32, loglevel: i32) -> i32;
fn trans_viscosity(n: i32) -> f64;
fn trans_thermalConductivity(n: i32) -> f64;
fn trans_getThermalDiffCoeffs(n: i32, ldt: i32, dt: *mut f64) -> i32; // Mass-averaged
}

use combustion::*;

pub fn check(model: &model::Model, state: &State) {
	let (species_names, species) = Species::new(&model.species);
	let transport_polynomials = species.transport_polynomials();

	let len = species.len();
	let file = std::ffi::CStr::from_bytes_with_nul(b"gri30.yaml\0").unwrap().as_ptr();
	let name = std::ffi::CStr::from_bytes_with_nul(b"gri30\0").unwrap().as_ptr();
	let phase = unsafe{thermo_newFromFile(file, name)};
	let cantera_species_names = (0..len).map(|k| {
		let mut specie = [0; 8];
		unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
		unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
	}).collect::<Box<_>>();
	assert!(unsafe{thermo_nSpecies(phase)} == len);
	assert!(state.amounts.len() == len && !state.amounts.iter().any(|&n| n<0.));
	let cantera_order = |o: &[f64]| (0..o.len()).map(|i| o[species_names.iter().position(|&s| s==cantera_species_names[i]).unwrap()]).collect::<Box<_>>();
	unsafe{thermo_setMoleFractions(phase, state.amounts.len(), cantera_order(&state.amounts).as_ptr(), 1)}; // /!\ Needs to be set before pressure
	dbg!(state.pressure_R * (K*NA), state.temperature);
	unsafe{thermo_setTemperature(phase, state.temperature)};
	unsafe{thermo_setPressure(phase, state.pressure_R * (K*NA))}; // /!\ Needs to be set after mole fractions
	let transport = unsafe{trans_newDefault(phase, 5)};
	let cantera_viscosity = unsafe{trans_viscosity(transport)};
	let cantera_thermal_conductivity  = unsafe{trans_thermalConductivity(transport)};
	//let ref binary_thermal_diffusion_coefficients = {let mut array = vec![0.; len*len]; unsafe{trans_getBinDiffCoeffs(transport, len as i32, array.as_mut_ptr())}; array};
	let ref cantera_mixture_mass_averaged_thermal_diffusion_coefficients = {
		let mut array = vec![0.; len];
		let iok = unsafe{trans_getThermalDiffCoeffs(transport, len as i32, array.as_mut_ptr())}; // /!\ Mass-averaged
		println!("{:?}", &array);
		assert!(iok == 0, "{}", iok);
		if false {
			array.into_boxed_slice()
		} else {
			if true { // Molar, N2
				"5.88830303e-04 9.74876528e-04 2.54222365e-04 1.48570420e-04 2.49489253e-04 2.11459327e-04 1.65222080e-04 1.64120812e-04 2.38722194e-04 2.71708310e-04 1.85500267e-04 1.85500267e-04 1.81398673e-04 1.97357482e-04 1.63943507e-04 1.32876379e-04 1.41676796e-04 1.40531852e-04 1.38019881e-04 1.38019881e-04 1.37553058e-04 1.39820441e-04 1.38417244e-04 1.37105272e-04 1.36668970e-04 1.25653611e-04 1.24639892e-04 2.05464735e-04 1.19597783e-04 1.19597783e-04 2.26754242e-04 2.67747169e-04 2.62329021e-04 1.99451789e-04 1.60259336e-04 1.62640124e-04 1.44415354e-04 1.31154960e-04 1.64801742e-04 1.61420566e-04 1.40439564e-04 1.39179409e-04 2.05460769e-04 1.31725863e-04 1.31725863e-04 1.31725863e-04 1.32333100e-04 1.64317769e-04 1.61220023e-04 9.72742481e-05 9.68451144e-05 1.19049476e-04 1.18523749e-04"
			} else { // Molar, Ar
				"6.01596776e-04 1.00226482e-03 2.55056423e-04 1.64626535e-04 2.50065057e-04 2.11551373e-04 1.62992833e-04 1.61816274e-04 2.39813147e-04 2.73446037e-04 1.85330651e-04 1.85330651e-04 1.81041369e-04 2.00619281e-04 1.62086393e-04 1.30003634e-04 1.39490186e-04 1.38275792e-04 1.35835729e-04 1.35835729e-04 1.35175678e-04 1.38074477e-04 1.36591137e-04 1.35202900e-04 1.34660395e-04 1.23590104e-04 1.22515786e-04 2.02913487e-04 1.16797862e-04 1.16797862e-04 2.27286765e-04 2.69056404e-04 2.63346607e-04 1.99152991e-04 1.58283213e-04 1.60618127e-04 1.41421583e-04 1.28295471e-04 1.62720202e-04 1.59699277e-04 1.38351491e-04 1.37018153e-04 2.02909207e-04 1.28907824e-04 1.28907824e-04 1.28907824e-04 1.29558757e-04 1.63492991e-04 1.36504022e-04 9.47262305e-05 9.42679458e-05 1.16211421e-04 1.15648787e-04"
			}
			.split(' ').map(|s| s.parse().unwrap()).collect::<Box<[f64]>>()
		}
	};
	let order = |o: &[f64]| (0..o.len()).map(|i| o[cantera_species_names.iter().position(|s| s==species_names[i]).unwrap()]).collect::<Box<_>>();
	let cantera_mixture_mass_averaged_thermal_diffusion_coefficients = order(&cantera_mixture_mass_averaged_thermal_diffusion_coefficients);

	let ref transport = combustion::transport::transport(&species.molar_mass, &transport_polynomials, state);

	use iter::{dot, zip, into::IntoCopied};
	let State{temperature, pressure_R, volume, amounts, ..} = state;
	let T = *temperature;
	let amount = pressure_R * volume / T;
	assert!(num::abs(amounts.iter().sum::<f64>()-amount) < 2e-15, "{:e}", num::abs(amounts.iter().sum::<f64>()-amount));
	let mole_fractions = iter::map(amounts.iter(), |n| n/amount);
	println!("Viscosity: {:.3}%", num::relative_error(transport.viscosity, cantera_viscosity)*100.);
	assert!(num::relative_error(transport.viscosity, cantera_viscosity) < 2e-5, "{:e}", num::relative_error(transport.viscosity, cantera_viscosity));

	println!("Thermal conductivity: {:.0}%", num::relative_error(transport.thermal_conductivity, cantera_thermal_conductivity)*100.);
	assert!(num::relative_error(transport.thermal_conductivity, cantera_thermal_conductivity) < 5e-2,
		"{:e}", num::relative_error(transport.thermal_conductivity, cantera_thermal_conductivity));
	/*let direct_thermal_conductivity = 1./2. * (
		dot(zip(mole_fractions.copied(), |k| species.thermal_conductivity(k, T))) +
		1. / dot(zip(mole_fractions.copied(), |k| 1. / species.thermal_conductivity(k, T)))
	);
	//println!("{:e}", num::relative_error(cantera_thermal_conductivity, direct_thermal_conductivity));
	//println!("{:e}", num::relative_error(transport.thermal_conductivity, direct_thermal_conductivity));
	assert!(num::relative_error(transport.thermal_conductivity, direct_thermal_conductivity) < 7e-4,
		"{:e}", num::relative_error(transport.thermal_conductivity, direct_thermal_conductivity));*/

	let molar_mass = dot(mole_fractions.iter().copied().zip(species.molar_mass.iter().copied()));
	let transport_mixture_mass_averaged_thermal_diffusion_coefficients = iter::eval(mole_fractions.len(), |k|
		(1. - mole_fractions[k]*species.molar_mass[k]/molar_mass) /
		dot(zip(mole_fractions.copied(), |j| if j != k { (pressure_R*K*NA) / transport_polynomials.binary_thermal_diffusion_coefficient(k, j, T) } else { 0. }))
	);

	for (specie, (&a, &b)) in		species_names.iter().zip(cantera_mixture_mass_averaged_thermal_diffusion_coefficients.iter().zip(transport_mixture_mass_averaged_thermal_diffusion_coefficients.iter())) {
		assert!(num::abs(a-b) < 1e-3 /*|| num::relative_error(a, b) < 0.*/, "{}: {:e} {:e}% {} {}", specie, num::abs(a-b), num::relative_error(a, b)*100., a, b);
	}
	/*use itertools::Itertools;
	println!("Diffusion coefficients: {}",
		cantera_mixture_mass_averaged_thermal_diffusion_coefficients.iter().zip(transport_mixture_mass_averaged_thermal_diffusion_coefficients.iter()).map(|(a,b)| format!("{:e} {:e}", a, b)).format("\n"));*/
	let e = cantera_mixture_mass_averaged_thermal_diffusion_coefficients.iter().zip(transport_mixture_mass_averaged_thermal_diffusion_coefficients.iter())
		.map(|(a,b)| num::relative_error(*a, *b)).filter(|e| e.is_finite()).reduce(f64::max).unwrap();
	println!("Diffusion coefficients: {:.0}%", e*100.);
	assert!(e < 0.2);
}
