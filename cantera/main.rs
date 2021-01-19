#![allow(mixed_script_confusables, non_snake_case, incomplete_features)]#![feature(type_ascription, array_map, non_ascii_idents, const_generics, const_evaluatable_checked)]
use combustion::*;
#[fehler::throws(Box<dyn std::error::Error>)] fn main() {
	let system = std::fs::read("CH4+O2.ron")?;
	test_transport_cantera(&Simulation::<35>::new(&system)?);
}

#[allow(dead_code)] fn test_transport_cantera<const S: usize>(Simulation{system, state, pressure_R, species_names, ..}: &Simulation<S>) where [(); S-1]: {
	let transport = system.transport(*pressure_R, &state);
	let cantera = {
		use itertools::Itertools;
		let mole_proportions = format!("{}", species_names.iter().zip(&state.amounts).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
		let mole_proportions = std::ffi::CString::new(mole_proportions).unwrap();
		use std::ptr::null;
		let ([mut viscosity, mut thermal_conductivity], mut species_len, mut specie_names, mut mixture_averaged_thermal_diffusion_coefficients) = ([0.; 2], 0, null(), null());
		unsafe {
			let pressure = pressure_R * (combustion::kB*combustion::NA);
			extern "C" { fn cantera(pressure: f64, temperature: f64, mole_proportions: *const std::os::raw::c_char,
				viscosity: &mut f64, thermal_conductivity: &mut f64,
				species_len: &mut usize, species: &mut *const *const std::os::raw::c_char, mixture_averaged_thermal_diffusion_coefficients: &mut *const f64); }
			cantera(pressure, state.temperature, mole_proportions.as_ptr(), &mut viscosity, &mut thermal_conductivity, &mut species_len, &mut specie_names, &mut mixture_averaged_thermal_diffusion_coefficients);
			let specie_names = iter::box_collect(std::slice::from_raw_parts(specie_names, species_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
			let order = |o:&[_]| iter::vec::eval(species_names, |s| o[specie_names.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}", s, species_names))]);
			let mixture_averaged_thermal_diffusion_coefficients = order(std::slice::from_raw_parts(mixture_averaged_thermal_diffusion_coefficients, species_len));
			Transport{viscosity, thermal_conductivity, mixture_averaged_thermal_diffusion_coefficients}
		}
	};
	dbg!(&transport, &cantera);
	dbg!((transport.viscosity-cantera.viscosity)/cantera.viscosity, (transport.thermal_conductivity-cantera.thermal_conductivity)/cantera.thermal_conductivity);
	assert!(f64::abs(transport.viscosity-cantera.viscosity)/cantera.viscosity < 0.03, "{}", transport.viscosity/cantera.viscosity);
	assert!(f64::abs(transport.thermal_conductivity-cantera.thermal_conductivity)/cantera.thermal_conductivity < 0.05, "{:?}", (transport.thermal_conductivity, cantera.thermal_conductivity));
}
