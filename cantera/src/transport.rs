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
	let ref mixture_averaged_thermal_diffusion_coefficients = {let mut array = vec![0.; len]; unsafe{trans_getThermalDiffCoeffs(transport, len as i32, array.as_mut_ptr())}; array};
	let ref transport = combustion::transport::transport(&species.molar_mass, &transport_polynomials, state);
	dbg!(&transport);
	dbg!(viscosity, transport.viscosity, num::relative_error(transport.viscosity, viscosity));
	dbg!(thermal_conductivity, transport.thermal_conductivity, num::relative_error(transport.thermal_conductivity, thermal_conductivity));
	if false {
		dbg!(mixture_averaged_thermal_diffusion_coefficients);
		dbg!(&transport.mixture_averaged_thermal_diffusion_coefficients);
		for (&a, &b) in mixture_averaged_thermal_diffusion_coefficients.iter().zip(transport.mixture_averaged_thermal_diffusion_coefficients.iter()) {
			dbg!(a, b, num::relative_error(a,b));
		}
	}
	assert!(num::relative_error(transport.viscosity, viscosity) < 0.03);
	assert!(num::relative_error(transport.thermal_conductivity, thermal_conductivity) < 0.05);
}
