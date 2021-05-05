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
//fn trans_getBinDiffCoeffs(n: i32, ld: i32, d: *mut f64) -> i32;
fn trans_getMixDiffCoeffs(n: i32, ldt: i32, dt: *mut f64) -> i32;
}

use combustion::*;

pub fn check(model: &model::Model, state: &State) {
	let pressure = state.pressure_R * (K*NA);
	let temperature = state.temperature;
	let (species_names, ref species) = Species::new(&model.species);

	let len = species.len();
	let file = std::ffi::CStr::from_bytes_with_nul(b"gri30.yaml\0").unwrap().as_ptr();
	let name = std::ffi::CStr::from_bytes_with_nul(b"gri30\0").unwrap().as_ptr();
	let phase = unsafe{thermo_newFromFile(file, name)};
	let _species_names = iter::eval(unsafe{thermo_nSpecies(phase)}, |k| {
		let mut specie = [0; 8];
		unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
		unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
	});
	let position = |i| _species_names.iter().position(|s| s==species_names[i]).unwrap();

	/*use std::io::BufRead;
	for line in std::fs::read("tdo").unwrap().lines() {
		let line = line.unwrap();
		use std::convert::TryInto;
		let [i,j, T⃰, δ⃰, Ω⃰22, A⃰, Ω⃰11]: [&str; 7] = (&*iter::box_collect(line.split('\t'))).try_into().unwrap();
		let [i,j] = [i,j].map(|i| position(i.parse().unwrap()));
		let check = |id, a, b| assert!(num::relative_error(a,b) < 2e-5, "{} {:e} {} {} {} {}", id, num::relative_error(a,b), i,j, a,b);
		check("T", species.T⃰(i,j, temperature), T⃰.parse().unwrap());
		check("δ", species.reduced_dipole_moment(i,j), δ⃰.parse().unwrap());
		check("Ω22", species.Ω⃰22(i,j, temperature), Ω⃰22.parse().unwrap());
		check("A", species.collision_integral(&transport::A⃰, i, j, temperature), A⃰.parse().unwrap());
		check("Ω11", species.Ω⃰11(i,j, temperature), Ω⃰11.parse().unwrap());
	}*/

	let ref transport_polynomials = species.transport_polynomials();
	let transport::Transport{viscosity, thermal_conductivity, mixture_diffusion_coefficients} =
		transport::transport(&species.molar_mass, &transport_polynomials, state);
	let mixture_diffusion_coefficients = iter::map(&*mixture_diffusion_coefficients, |cP| cP / pressure);

	let _order = |o: &[f64]| iter::eval(o.len(), |i| o[species_names.iter().position(|&s| s==_species_names[i]).unwrap()]);
	unsafe{thermo_setMoleFractions(phase, state.amounts.len(), _order(&state.amounts).as_ptr(), 1)}; // /!\ Needs to be set before pressure
	unsafe{thermo_setTemperature(phase, temperature)};
	unsafe{thermo_setPressure(phase, pressure)}; // /!\ Needs to be set after mole fractions
	let transport = unsafe{trans_newDefault(phase, 5)};
	let _viscosity = unsafe{trans_viscosity(transport)};
	let _thermal_conductivity  = unsafe{trans_thermalConductivity(transport)};
	//let ref _binary_thermal_diffusion_coefficients = {let mut array = vec![0.; len*len]; unsafe{trans_getBinDiffCoeffs(transport, len as i32, array.as_mut_ptr())}; array};
	//let ref _binary_thermal_diffusion_coefficients = iter::map(_binary_thermal_diffusion_coefficients, |c_P| pressure*c_P);
	let ref _mixture_diffusion_coefficients = {
		let mut array = vec![0.; len];
		assert!(unsafe{trans_getMixDiffCoeffs(transport, len as i32, array.as_mut_ptr())} == 0); // Mass-averaged
		if true {
			array.into_boxed_slice()
		} else { // N2
			iter::box_collect("5.88830303e-04 9.74876528e-04 2.54222365e-04 1.48570420e-04 2.49489253e-04 2.11459327e-04 1.65222080e-04 1.64120812e-04 2.38722194e-04 2.71708310e-04 1.85500267e-04 1.85500267e-04 1.81398673e-04 1.97357482e-04 1.63943507e-04 1.32876379e-04 1.41676796e-04 1.40531852e-04 1.38019881e-04 1.38019881e-04 1.37553058e-04 1.39820441e-04 1.38417244e-04 1.37105272e-04 1.36668970e-04 1.25653611e-04 1.24639892e-04 2.05464735e-04 1.19597783e-04 1.19597783e-04 2.26754242e-04 2.67747169e-04 2.62329021e-04 1.99451789e-04 1.60259336e-04 1.62640124e-04 1.44415354e-04 1.31154960e-04 1.64801742e-04 1.61420566e-04 1.40439564e-04 1.39179409e-04 2.05460769e-04 1.31725863e-04 1.31725863e-04 1.31725863e-04 1.32333100e-04 1.64317769e-04 1.61220023e-04 9.72742481e-05 9.68451144e-05 1.19049476e-04 1.18523749e-04".split(' ').map(|s| s.parse().unwrap()))
		}
	};
	let order = |o: &[f64]| iter::eval(o.len(), |i| o[position(i)]);
	//let order2 = |o: &[f64]| iter::box_collect((0..len).map(|i| (0..len).map(move |j| o[position(i)*len+position(j)])).flatten());
	//let _binary_thermal_diffusion_coefficients = order2(&_binary_thermal_diffusion_coefficients);
	let _mixture_diffusion_coefficients = order(&_mixture_diffusion_coefficients);

	println!("Viscosity: {:.3}%", num::relative_error(viscosity, _viscosity)*100.);
	println!("Thermal conductivity: {:.0}%", num::relative_error(thermal_conductivity, _thermal_conductivity)*100.);
	assert!(num::relative_error(viscosity, _viscosity) < 2e-5, "{:e}", num::relative_error(viscosity, _viscosity));
	assert!(num::relative_error(thermal_conductivity, _thermal_conductivity) < 6e-2, "{:e}", num::relative_error(thermal_conductivity, _thermal_conductivity));

	//let binary_thermal_diffusion_coefficients = iter::box_collect((0..len).map(|i| (0..len).map(move |j| species.binary_thermal_diffusion_coefficient(i,j, temperature))).flatten());
	/*let binary_thermal_diffusion_coefficients = iter::box_collect((0..len).map(|i| (0..len).map(move |j| transport_polynomials.binary_thermal_diffusion_coefficient(i,j, temperature))).flatten());
	for (specie, (&a, &b)) in	species_names.iter().zip(binary_thermal_diffusion_coefficients.iter().zip(_binary_thermal_diffusion_coefficients.iter())) {
		assert!(num::relative_error(a, b) < 4e-3, "{}: {} {} {:e} {:e}", specie, a, b, num::abs(a-b), num::relative_error(a, b));
	}*/

	/*let amount = pressure_R * volume / temperature;
	let mole_fractions = iter::map(amounts.iter(), |n| n/amount);
	use iter::{dot, zip, into::IntoCopied};
	let molar_mass = dot(mole_fractions.iter().copied().zip(species.molar_mass.iter().copied()));
	let mixture_mass_averaged_thermal_diffusion_coefficients = iter::eval(mole_fractions.len(), |k|
		(1. - mole_fractions[k]*species.molar_mass[k]/molar_mass) /
		dot(zip(mole_fractions.copied(), |j| if j != k { /*(pressure_R*K*NA)*/1. / species.binary_thermal_diffusion_coefficient(k, j, T) } else { 0. }))
	);*/

	for (_specie, (&a, &b)) in	species_names.iter().zip(mixture_diffusion_coefficients.iter().zip(_mixture_diffusion_coefficients.iter())) {
		assert!(num::relative_error(a, b) < 1e-2, "{:e}", num::relative_error(a, b));
	}
	/*use itertools::Itertools;
	println!("Diffusion coefficients: {}",
		cantera_mixture_mass_averaged_thermal_diffusion_coefficients.iter().zip(transport_mixture_mass_averaged_thermal_diffusion_coefficients.iter()).map(|(a,b)| format!("{:.2e} {:.2e}", a, b)).format("\n"));
	let e = cantera_mixture_mass_averaged_thermal_diffusion_coefficients.iter().zip(transport_mixture_mass_averaged_thermal_diffusion_coefficients.iter())
		.map(|(a,b)| num::relative_error(*a, *b)).filter(|e| e.is_finite()).reduce(f64::max).unwrap();
	println!("Diffusion coefficients: {:.0}%", e*100.);
	assert!(e < 0.2);
	println!("{}", transport_mixture_mass_averaged_thermal_diffusion_coefficients.iter()).map(|(a,b)| format!("{:.2e} {:.2e}", a, b)).format("\n"));*/
}
