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
fn trans_getMixDiffCoeffs(n: i32, ldt: i32, dt: *mut f64) -> i32;
}

pub fn check(file: &str) {
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&file).unwrap()).unwrap()).unwrap();
	let model = yaml::parse(&model).unwrap();
	let state = combustion::initial_state(&model);
	let (species_names, ref species) = combustion::Species::new(&model.species);
	let pressure = 101325.;
	let temperature = 1000.;

	let len = species.len();
	let phase = unsafe{let file = std::ffi::CString::new(file).unwrap(); thermo_newFromFile(file.as_ptr(), std::ffi::CStr::from_bytes_with_nul(b"\0").unwrap().as_ptr())};
	let _species_names = iter::eval(unsafe{thermo_nSpecies(phase)}, |k| {
		let mut specie = [0; 8];
		unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
		unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
	});
	let position = |i| _species_names.iter().position(|s| s==species_names[i]).unwrap();
	let ref transport_polynomials = species.transport_polynomials();
	let transport::Transport{viscosity, thermal_conductivity, mixture_diffusion_coefficients} = transport::transport(&species.molar_mass, &transport_polynomials, &state);
	let mixture_diffusion_coefficients = iter::map(&*mixture_diffusion_coefficients, |cP| cP / pressure);
	let _order = |o: &[f64]| iter::eval(o.len(), |i| o[species_names.iter().position(|&s| s==_species_names[i]).unwrap()]);
	unsafe{thermo_setMoleFractions(phase, state.amounts.len(), _order(&state.amounts).as_ptr(), 1)}; // /!\ Needs to be set before pressure
	unsafe{thermo_setTemperature(phase, temperature)};
	unsafe{thermo_setPressure(phase, pressure)}; // /!\ Needs to be set after mole fractions
	let transport = unsafe{trans_newDefault(phase, 0)};
	let _viscosity = unsafe{trans_viscosity(transport)};
	let _thermal_conductivity  = unsafe{trans_thermalConductivity(transport)};
	let ref _mixture_diffusion_coefficients = {let mut array = vec![0.; len]; assert!(unsafe{trans_getMixDiffCoeffs(transport, len as i32, array.as_mut_ptr())} == 0); array};
	let order = |o: &[f64]| iter::eval(o.len(), |i| o[position(i)]);
	let _mixture_diffusion_coefficients = order(&_mixture_diffusion_coefficients);
	println!("Viscosity: {:.0e}", num::relative_error(viscosity, _viscosity));
	println!("Thermal conductivity: {:.0e}", num::relative_error(thermal_conductivity, _thermal_conductivity));
	let e = mixture_diffusion_coefficients.iter().zip(_mixture_diffusion_coefficients.iter()).map(|(a,b)| num::relative_error(*a, *b)).filter(|e| e.is_finite()).reduce(f64::max).unwrap();
	println!("Mixture diffusion coefficients: {:.0e}", e);
	assert!(num::relative_error(viscosity, _viscosity) < 5e-4, "{:e}", num::relative_error(viscosity, _viscosity));
	assert!(num::relative_error(thermal_conductivity, _thermal_conductivity) < 0.06, "{:e}", num::relative_error(thermal_conductivity, _thermal_conductivity));
	assert!(e < 0.01);
	use itertools::Itertools;
	println!("{}", mixture_diffusion_coefficients.iter().format(" "));
}
