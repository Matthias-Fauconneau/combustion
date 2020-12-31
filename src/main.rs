#![feature(type_ascription, array_map, non_ascii_idents)]#![allow(mixed_script_confusables,non_snake_case)]
extern "C" {
fn cantera(pressure: f64, temperature: f64, mole_proportions: *const std::os::raw::c_char,
									viscosity: &mut f64, thermal_conductivity: &mut f64,
									species_len: &mut usize, species: &mut *const *const std::os::raw::c_char, multicomponent_diffusion_coefficients: &mut *const f64);

}

#[fehler::throws(Box<dyn std::error::Error>)] fn main() {
	let system = std::fs::read("CH4+O2.ron")?;
	const S : usize = 35; // Number of species
	type Simulation<'t> = combustion::Simulation::<'t, S>;
	let Simulation{system, state: combustion::State{temperature, amounts}, pressure_r, species, ..} = Simulation::new(&system)?;
	let pressure = pressure_r * (combustion::kB*combustion::NA);
	let ([viscosity, _thermal_conductivity], ref _thermal_diffusion_coefficients) = {
		use itertools::Itertools;
		let mole_proportions = format!("{}", species.iter().zip(&amounts).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
		let mole_proportions = std::ffi::CString::new(mole_proportions).unwrap();
		use std::ptr::null;
		let ([mut viscosity, mut thermal_conductivity], mut species_len, mut specie_names, mut thermal_diffusion_coefficients) = ([0.; 2], 0, null(), null());//:([f64; 2], usize, *const *const std::os::raw::c_char, *const f64)
		unsafe {
			dbg!(pressure, temperature, &mole_proportions);
			cantera(pressure, temperature, mole_proportions.as_ptr(), &mut viscosity, &mut thermal_conductivity, &mut species_len, &mut specie_names, &mut thermal_diffusion_coefficients);
			let specie_names = iter::box_collect(std::slice::from_raw_parts(specie_names, species_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
			let order = |o:&[_]| iter::vec::eval(species, |s| o[specie_names.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}",s,species))]);
			let thermal_diffusion_coefficients = std::slice::from_raw_parts(thermal_diffusion_coefficients, species_len);
			([viscosity, thermal_conductivity], order(thermal_diffusion_coefficients))
		}
	};
	//dbg!(([viscosity, thermal_conductivity], /*thermal_diffusion_coefficients*/));
	//dbg!(system.transport(pressure, temperature, amounts));
	assert_eq!(system.transport(pressure, temperature, &amounts), viscosity);
}
