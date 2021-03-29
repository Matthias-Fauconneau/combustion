//fn dot(iter: impl IntoIterator<Item=(f64, f64)>) -> f64 { iter.into_iter().map(|(a,b)| a*b).sum() }
use super::*;

pub fn check(model: model::Model, Simulation{state, ..}: &Simulation) {
	let (_species_names, species) = Species::new(model.species);
	let transport_polynomials = species.transport_polynomials();

	let len = species.len();
	let file = std::ffi::CStr::from_bytes_with_nul(b"gri30.yaml\0").unwrap().as_ptr();
	let name = std::ffi::CStr::from_bytes_with_nul(b"gri30\0").unwrap().as_ptr();
	let phase = unsafe{thermo_newFromFile(file, name)};
	assert!(unsafe{thermo_nSpecies(phase)} == len);
	assert!(state.amounts.len() == len && !state.amounts.iter().any(|&n| n<0.));
	unsafe{thermo_setMoleFractions(phase, state.amounts.len(), state.amounts.as_ptr(), 1)}; // /!\ Needs to be set before pressure
	dbg!(state.pressure * NA, state.temperature);
	unsafe{thermo_setTemperature(phase, state.temperature)};
	unsafe{thermo_setPressure(phase, state.pressure * NA)}; // /!\ Needs to be set after mole fractions
	let transport = unsafe{trans_newDefault(phase, 5)};
	let viscosity = unsafe{trans_viscosity(transport)};
	let thermal_conductivity  = unsafe{trans_thermalConductivity(transport)};
	let mixture_averaged_thermal_diffusion_coefficients = {let mut array = vec![0.; len]; unsafe{trans_getThermalDiffCoeffs(transport, len as i32, array.as_mut_ptr())}; array};
	let transport = combustion::transport::transport(&species.molar_mass, &transport_polynomials, state);
	dbg!(&transport, viscosity, thermal_conductivity, mixture_averaged_thermal_diffusion_coefficients);
	dbg!((transport.viscosity-viscosity)/viscosity, (transport.thermal_conductivity-thermal_conductivity)/thermal_conductivity);
	assert!(f64::abs(transport.viscosity-viscosity)/viscosity < 0.03, "{}", transport.viscosity/viscosity);
	assert!(f64::abs(transport.thermal_conductivity-thermal_conductivity)/thermal_conductivity < 0.05, "{:?}", (transport.thermal_conductivity, thermal_conductivity));
}
