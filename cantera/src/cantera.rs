//fn dot(iter: impl IntoIterator<Item=(f64, f64)>) -> f64 { iter.into_iter().map(|(a,b)| a*b).sum() }
use {itertools::Itertools, super::*};

pub fn check(model: Model, Simulation{species_names, state, ..}: &Simulation) {
	let len = model.len();
	let file = std::ffi::CStr::from_bytes_with_nul(b"gri30.yaml\0").unwrap().as_ptr();
	let name = std::ffi::CStr::from_bytes_with_nul(b"gri30\0").unwrap().as_ptr();
	let phase = unsafe{thermo_newFromFile(file, name)};
	assert!(unsafe{thermo_nSpecies(phase)} == len);
	let kinetics = unsafe{kin_newFromFile(file, name, phase, 0, 0, 0, 0)};
	assert!(state.amounts.len() == len && !state.amounts.iter().any(|&n| n<0.));
	unsafe{thermo_setMoleFractions(phase, state.amounts.len(), state.amounts.as_ptr(), 1)}; // /!\ Needs to be set before pressure
	dbg!(state.pressure * NA, state.temperature);
	unsafe{thermo_setTemperature(phase, state.temperature)};
	unsafe{thermo_setPressure(phase, state.pressure * NA)}; // /!\ Needs to be set after mole fractions
	let mut specie_production_rates = vec![0.; len];
	unsafe{kin_getNetProductionRates(kinetics, len, specie_production_rates.as_mut_ptr())};
	let specie_production_rates = specie_production_rates.iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>(); // kmol -> mol
	/*let mut Cp_R = vec![0.; len];
	unsafe{thermo_getCp_R(phase, len, Cp_R.as_mut_ptr())};
	let volume = state.volume;
	let ref concentrations = state.amounts.iter().map(|n| n / volume).collect::<Box<_>>();
	let ΣCCp = dot(concentrations.iter().copied().zip(Cp_R.iter().copied()));
	let mut H_RT = vec![0.; len];
	unsafe{thermo_getEnthalpies_RT(phase, len, H_RT.as_mut_ptr())};
	let dtT_T = -1./ΣCCp * dot(specie_production_rates.iter().copied().zip(H_RT.iter().copied()));
	let T = state.temperature;
	let dtT = dtT_T*T;*/
	println!("{}", specie_production_rates.iter().zip(species_names.iter()).filter(|(&v,_)| v!=0.).map(|(v,name)| format!("{}: {:.6e}",name, v)).format(", "));
}
