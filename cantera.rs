#![feature(format_args_capture,iter_is_partitioned)]#![allow(non_snake_case)]
use std::os::raw::c_char;
#[link(name = "cantera")]
extern "C" {
fn thermo_newFromFile(file_name: *const c_char, phase_name: *const c_char) -> i32;
fn thermo_nSpecies(n: i32) -> usize;
fn thermo_setTemperature(n: i32, t: f64) -> i32;
fn thermo_setMoleFractions(n: i32, len: usize, x: *const f64, norm: i32) -> i32;
fn thermo_getSpeciesName(n: i32, m: usize, len: usize, buffer: *mut c_char) -> i32;
fn thermo_setPressure(n: i32, p: f64) -> i32;
fn kin_newFromFile(file_name: *const c_char, phase_name: *const c_char, reactingPhase: i32, neighbor0: i32, neighbor1: i32, neighbor2: i32, neighbor3: i32) -> i32;
fn kin_getNetProductionRates(n: i32, len: usize, w_dot: *mut f64) -> i32;
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let path = std::env::args().skip(1).next().unwrap();
	let model = yaml_model::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path)?)?)?;
	let model = yaml_model::parse(&model);
	use chemistry::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	let reactions = map(&*model.reactions, |r| Reaction::new(species_names, r));
	let active = {
		let active = map(0..species.len()-1, |k| reactions.iter().any(|Reaction{net,..}| net[k] != 0));
		assert!(active.iter().is_partitioned(|&active| active),
								"Species must be partionned with all active species first and all inert species last so that the state only consists of active species when solving kinetics: {active:?}");
		active.iter().position(|active| !active).unwrap_or(species.len()-1)
	};

	let file = std::ffi::CString::new(path).unwrap();
	let phase_name_cstr_ptr = std::ffi::CStr::from_bytes_with_nul(b"gri30\0").unwrap().as_ptr();
	let phase = unsafe{thermo_newFromFile(file.as_c_str().as_ptr(), phase_name_cstr_ptr)};
	let cantera_species_names = map(0..species.len(), |k| {
		let mut specie = [0; 8];
		unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
		unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
	});
	assert!(unsafe{thermo_nSpecies(phase)} == species.len());
	let kinetics = unsafe{kin_newFromFile(file.as_c_str().as_ptr(), phase_name_cstr_ptr, phase, 0, 0, 0, 0)};
	let cantera_order = |o: &[f64]| (0..o.len()).map(|i| o[species_names.iter().position(|&s| s==cantera_species_names[i]).unwrap()]).collect::<Box<_>>();

	let ref state = initial_state(&model);
	use {iter::map, ast::{wrap, let_}};
	let rates = wrap(reaction::rates(species.len(), &species.thermodynamics[0..active], &reactions));
	assert!(state.volume == 1.);
	let test = |state: &State| {
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let active_amounts = &amounts[0..amounts.len()-1];
		let_!{ [_energy_rate_RT, rates @ ..] = &*rates(&[&[*pressure_R, total_amount, *T], active_amounts].concat()) => {
		//println!("{}, HRR: {:.3e}", rates.iter().zip(&**species_names).format_with(", ", |(rate, name), f| f(&format!("{name}: {rate:.0}").to_string())), NA * kB * T * -energy_rate_RT);

		unsafe{thermo_setMoleFractions(phase, state.amounts.len(), cantera_order(&state.amounts).as_ptr(), 1)}; // /!\ Needs to be set before pressure
		unsafe{thermo_setTemperature(phase, state.temperature)};
		unsafe{thermo_setPressure(phase, state.pressure_R * (kB*NA))}; // /!\ Needs to be set after mole fractions
		let order = |o: &[f64]| map(&**species_names, |specie| o[cantera_species_names.iter().position(|s| s==specie).unwrap()]);
		let cantera_rates = {
			let mut rates = vec![0.; species.len()];
			unsafe{kin_getNetProductionRates(kinetics, rates.len(), rates.as_mut_ptr())};
			order(&map(rates, |c| c*1000.)) // kmol -> mol
		};
		rates.iter().zip(&cantera_rates[0..species.len()-1]).map(|(&a,&b)| num::relative_error(a,b)).reduce(f64::max).unwrap()
	}}};
	let mut random= rand::thread_rng();
	let mut max = 0.;
	for i in 0.. {
		let volume = 1.;
		use rand::Rng;
		let temperature = random.gen_range(800. .. 1200.);
		//let temperature = random.gen_range(1000. .. 1001.);
		//let temperature = 1000.;
		let pressure_R = random.gen_range(0. .. 10e6)/(kB*NA);
		//let pressure_R = random.gen_range(101325. .. 101326.)/(kB*NA);
		//let pressure_R = 101325./(kB*NA);
		let amount = pressure_R * volume / temperature;
		let amount_proportions = map(0..species.len(), |_| random.gen());
		//let amount_proportions = map(0..species.len(), |_| random.gen_range(1./(species.len() as f64) .. 1./(species.len() as f64-1.)));
		//let amount_proportions = map(0..species.len(), |_| 1./(species.len() as f64));
		let amounts = map(&*amount_proportions, |amount_proportion| amount * amount_proportion/amount_proportions.iter().sum::<f64>());
		let e = test(&State{temperature, pressure_R, volume, amounts});
		if e > max { max = e; println!("{i} {e:.0e}"); }
	}
	//use itertools::Itertools; println!("{}", species_names.iter().zip(rates.iter().zip(&cantera_rates[0..species.len()-1])).format_with(", ", |(name, (&a,&b)), f| f(&format!("{name}: {:.0e}", num::relative_error(a,b)).to_string())));
	Ok(())
}
