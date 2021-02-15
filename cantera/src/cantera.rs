extern "C" {
pub fn reaction(pressure: &mut f64, temperature: &mut f64, mole_proportions: *const std::os::raw::c_char,
													species_len: &mut usize, species: &mut *const *const std::os::raw::c_char,
													rate_time: f64,
														net_productions_rates: &mut *const f64,
														reactions_len: &mut usize,
															equations: &mut *const *const std::os::raw::c_char,
															equilibrium_constants: &mut *const f64,
															forward: &mut *const f64,
															reverse: &mut *const f64,
													state_time: f64,
													concentrations: &mut *const f64);
}
extern "C" {
pub fn transport(pressure: f64, temperature: f64, mole_proportions: *const std::os::raw::c_char,
														viscosity: &mut f64, thermal_conductivity: &mut f64,
														species_len: &mut usize, species: &mut *const *const std::os::raw::c_char, mixture_averaged_thermal_diffusion_coefficients: &mut *const f64);
}
