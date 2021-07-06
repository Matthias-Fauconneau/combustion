#![feature(format_args_capture,iter_partition_in_place,array_map)]#![allow(non_snake_case,non_upper_case_globals)]
mod yaml; mod device;
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
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path)?)?)?;
	let model = yaml::parse(&model);
	use combustion::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	let reactions = map(&*model.reactions, |r| Reaction::new(species_names, r));
	let rates = reaction::rates(&species.thermodynamics, &reactions);
	let rates = device::assemble(&rates);

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

	let ref state = initial_state(&model);
	use {iter::map, ast::let_};
	assert!(state.volume == 1.);
	let test = move |state: &State| {
		use itertools::Itertools;

		let State{temperature, pressure_R, amounts, ..} = state;

		let total_amount = amounts.iter().sum::<f64>();
		let active_amounts = &amounts[0..amounts.len()-1];
		let input = iter::box_([total_amount, *temperature].iter().chain(active_amounts).copied().map(|input| vec![input as _; 1].into_boxed_slice()));
		let_!{ [_energy_rate_RT, rates @ ..] = &*rates(&[*pressure_R as _], &map(&*input, |input| &**input)).unwrap() => {
		let rates = map(&*rates, |v| v[0] as _);
		assert!(rates.len() == species.len()-1, "{}", rates.len());
		if false { println!("{}", rates.iter().format(" ")); }

		let cantera_order = |o: &[f64]| (0..o.len()).map(|i| o[species_names.iter().position(|&s| s==cantera_species_names[i]).unwrap()]).collect::<Box<_>>();
		unsafe{thermo_setMoleFractions(phase, amounts.len(), cantera_order(&amounts).as_ptr(), 1)}; // /!\ Needs to be set before pressure
		unsafe{thermo_setTemperature(phase, *temperature)};
		unsafe{thermo_setPressure(phase, pressure_R * (kB*NA))}; // /!\ Needs to be set after mole fractions
		let cantera = {
			let mut rates = vec![0.; species.len()];
			unsafe{kin_getNetProductionRates(kinetics, rates.len(), rates.as_mut_ptr())};
			let order = |o: &[f64]| map(&species_names[0..species.len()-1], |specie| o[cantera_species_names.iter().position(|s| s==specie).unwrap()]);
			order(&map(rates, |c| c*1000.)) // kmol -> mol
		};
		assert!(cantera.len() == species.len()-1);
		if false { println!("{}", cantera.iter().format(" ")); }

		let amounts = cantera_species_names.iter().map(|specie| amounts[species_names.iter().position(|s| s==specie).unwrap()]).format(" ");
		if true { rates.iter().zip(&*cantera).map(|(&a,&b)| num::relative_error(a,b)).reduce(f64::max).unwrap() }
		else {
			let nekrk = std::process::Command::new("../nekRK/build/main").current_dir("../nekRK").args(["Serial","1","1","0","gri30",&amounts.to_string(), &temperature.to_string(), &(pressure_R*(kB*NA)).to_string()]).output().unwrap();
			let nekrk = nekrk.stdout;
			let nekrk = std::str::from_utf8(&nekrk).unwrap();
			let nekrk : std::collections::HashMap<&str,f64> = nekrk.split(", ").map(|entry| {
				let entry:Box<_> = entry.split(": ").collect::<Box<_>>().try_into().unwrap();
				let [key, value]:[&str;2] = *entry;
				(key, value.trim().parse().expect(value))
			}).collect();
			let nekrk = map(&species_names[0..species.len()-1], |specie| nekrk[specie]); // Already in mol
			rates.iter().zip(&*cantera).zip(&*nekrk).map(|((&a,&b), &c)| f64::max(num::relative_error(a,b), num::relative_error(b,c))).reduce(f64::max).unwrap()
		}
	}}};
	let mut random= rand::thread_rng();
	let mut max = 0.;
	for i in 0.. {
		let volume = 1.;
		use rand::Rng;
		let temperature = random.gen_range(1000. .. 1200.);
		let pressure = random.gen_range(0. .. 1e5);
		let pressure_R = pressure/(kB*NA); //random.gen_range(0. .. 10e6)/(kB*NA);
		let amount = pressure_R * volume / temperature;
		let amount_proportions = map(0..species.len(), |_| random.gen());
		let amounts = map(&*amount_proportions, |amount_proportion| amount * amount_proportion/amount_proportions.iter().sum::<f64>());
		let e = test(&State{temperature, pressure_R, volume, amounts});
		//println!("{temperature} {pressure} {e}");
		assert!(e < 1e-2, "T {temperature} P {pressure} E {e} {max}");
		if e > max { max = e; println!("{i} {e:.0e}"); }
	}
	Ok(())
}
