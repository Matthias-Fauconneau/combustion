#![feature(type_ascription, array_map, non_ascii_idents)]#![allow(mixed_script_confusables, non_snake_case)]
extern "C" {
fn cantera(pressure: f64, temperature: f64, mole_proportions: *const std::os::raw::c_char,
									viscosity: &mut f64, thermal_conductivity: &mut f64,
									species_len: &mut usize, species: &mut *const *const std::os::raw::c_char, mixture_averaged_thermal_diffusion_coefficients: &mut *const f64);

}

#[fehler::throws(Box<dyn std::error::Error>)] fn main() {
	//trace::rstack_self()?; trace::sigfpe();
	let system = std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let Simulation{system, state: combustion::State{amounts, ..}, pressure_R, species, ..} = Simulation::<35>::new(&system)?;
	let temperature = 273.15+500.;
	let transport = system.transport(pressure_R, temperature, &amounts);
	//trace::unmask_SSE_exceptions();
	let cantera = {
		use itertools::Itertools;
		let mole_proportions = format!("{}", species.iter().zip(&amounts).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
		let mole_proportions = std::ffi::CString::new(mole_proportions).unwrap();
		use std::ptr::null;
		let ([mut viscosity, mut thermal_conductivity], mut species_len, mut specie_names, mut mixture_averaged_thermal_diffusion_coefficients) = ([0.; 2], 0, null(), null());
		unsafe {
			let pressure = pressure_R * (combustion::kB*combustion::NA);
			//dbg!(/*pressure,*/ temperature, &mole_proportions);
			cantera(pressure, temperature, mole_proportions.as_ptr(), &mut viscosity, &mut thermal_conductivity, &mut species_len, &mut specie_names, &mut mixture_averaged_thermal_diffusion_coefficients);
			let specie_names = iter::box_collect(std::slice::from_raw_parts(specie_names, species_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
			let order = |o:&[_]| iter::vec::eval(species, |s| o[specie_names.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}",s,species))]);
			let mixture_averaged_thermal_diffusion_coefficients = order(std::slice::from_raw_parts(mixture_averaged_thermal_diffusion_coefficients, species_len));
			Transport{viscosity, thermal_conductivity, mixture_averaged_thermal_diffusion_coefficients}
		}
	};
	dbg!(&transport, &cantera);
	dbg!((transport.viscosity-cantera.viscosity)/cantera.viscosity, (transport.thermal_conductivity-cantera.thermal_conductivity)/cantera.thermal_conductivity);
	assert!(f64::abs(transport.viscosity-cantera.viscosity)/cantera.viscosity < 0.04, "{}", transport.viscosity/cantera.viscosity);
	assert!(f64::abs(transport.thermal_conductivity-cantera.thermal_conductivity)/cantera.thermal_conductivity < 0.06, "{:?}", (transport.thermal_conductivity, cantera.thermal_conductivity));
}
