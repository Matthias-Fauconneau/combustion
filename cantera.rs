#![feature(format_args_capture,array_map)]#![allow(non_snake_case,non_upper_case_globals)]
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
fn kin_getNetRatesOfProgress(n: i32, len: usize, rates: *mut f64) -> i32;
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let mechanism_file_path = std::env::args().skip(1).next().unwrap();
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&mechanism_file_path)?)?)?;
	let model = yaml::parse(&model);
	use combustion::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	//eprintln!("{species_names:?}");
	let reactions = map(&*model.reactions, |r| Reaction::new(species_names, r));
	let rates = reaction::rates(&species.thermodynamics, &reactions);
	let rates = device::assemble(&rates);

	let file = std::ffi::CString::new(mechanism_file_path.as_str()).unwrap();
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
	use {iter::{list, map}, ast::let_};
	assert!(state.volume == 1.);
	let test = move |state: &State| {
		use itertools::Itertools;

		let State{temperature, pressure_R, amounts, ..} = state;

		let total_amount = amounts.iter().sum::<f64>();
		let active_amounts = &amounts[0..amounts.len()-1];
		let input = list([total_amount, *temperature].iter().chain(active_amounts).copied().map(|input| vec![input as _; 1].into_boxed_slice()));
		let_!{ [_energy_rate_RT, rates @ ..] = &*rates(&[*pressure_R as _], &map(&*input, |input| &**input)).unwrap() => {
		let rates = map(&*rates, |v| v[0] as _);
		assert!(rates.len() == species.len()-1, "{}", rates.len());
		if false { println!("{}", rates.iter().format(" ")); }

		let cantera_order = |o: &[f64]| (0..o.len()).map(|i| o[species_names.iter().position(|&s| s==cantera_species_names[i]).unwrap()]).collect::<Box<_>>();
		let amount_fractions = map(&**amounts, |n| n/total_amount);
		unsafe{thermo_setMoleFractions(phase, amount_fractions.len(), cantera_order(&amount_fractions).as_ptr(), 1)}; // /!\ Needs to be set before pressure
		unsafe{thermo_setTemperature(phase, *temperature)};
		unsafe{thermo_setPressure(phase, pressure_R * (kB*NA))}; // /!\ Needs to be set after mole fractions
		let cantera = if false { // Species
			let mut rates = vec![0.; species.len()];
			unsafe{kin_getNetProductionRates(kinetics, rates.len(), rates.as_mut_ptr())};
			let order = |o: &[f64]| map(&species_names[0..species.len()-1], |specie| o[cantera_species_names.iter().position(|s| s==specie).unwrap()]);
			let species_rates = order(&map(rates, |c| c*1000.)); // kmol -> mol
			assert!(species_rates.len() == species.len()-1);
			species_rates
		} else { // Reactions
			let mut rates = vec![0.; reactions.len()];
			unsafe{kin_getNetRatesOfProgress(kinetics, rates.len(), rates.as_mut_ptr())};
			map(rates, |c| c*1000.) // kmol -> mol
		};
		if true{ println!("{:.0}", cantera.iter().format(" ")); }

		let error = |a,b| {
			let m=f64::abs(f64::max(a,b));
			if m < 1e2 { return 0.; }
			let threshold=1e6; (if m < threshold { m/threshold } else { 1. }) * num::relative_error(a,b)
		};
		//let error = |a,b| num::relative_error(a,b);
		if false {
			let max = rates.iter().zip(&*cantera).map(|(&a,&b)| error(a,b)).enumerate().reduce(|a,b| if a.1 > b.1 { a } else { b }).unwrap();
			if max.1 > 1e-4 {
				println!("{}", rates.iter().format(" "));
				println!("{} {}", rates[max.0], cantera[max.0]);
			}
			max
		}
		else {
			//amount_fractions is in split active/inert indexing (for Rust and nekRK) (TODO: parse amounts[specie name] in C++ to prevent such errors)
			//let amount_fractions = cantera_species_names.iter().map(|specie| amount_fractions[species_names.iter().position(|s| s==specie).unwrap()]).format(" ").to_string(); // Cantera order
			//println!("cd ../nekRK && ../nekRK/build/main Serial 1 1 0 gri30 {} {} {}",temperature,pressure_R*(kB*NA),amount_fractions);
			//std::fs::write("/var/tmp/main.c", std::process::Command::new("../nekRK/main.py").arg(&mechanism_file_path).output().unwrap().stdout).unwrap(); // FIXME: don't do it everytime since OCCA takes forever to pass the code on to nvcc
			let nekrk = std::process::Command::new("../nekRK/build/main").current_dir("../nekRK")
				.args(["/var/tmp/main.c","Serial","1","1","0",&(pressure_R*(kB*NA)).to_string(),&temperature.to_string(),&amount_fractions.iter().format(" ").to_string()]).output().unwrap();
			assert!(nekrk.stderr.is_empty(), "{}", std::str::from_utf8(&nekrk.stderr).unwrap());
			let stdout = nekrk.stdout;
			let stdout = std::str::from_utf8(&stdout).unwrap();
			let lines = list(stdout.split("\n"));
			let ref nekrk_species_names = list(lines[0].trim().split(" "));
			assert!(species_names == nekrk_species_names, "\n{species_names:?}\n{nekrk_species_names:?}");
			let nekrk = map(lines[1].trim().split(" "), |r| r.trim().parse().unwrap());
			let max = cantera.iter().zip(&*nekrk).map(|(&a,&b)| error(a,b)).enumerate().reduce(|a,b| if a.1 > b.1 { a } else { b }).unwrap();
			if max.1 > 1e-4 {
				println!("{:.0}", nekrk.iter().format(" "));
				println!("{:>6}", cantera.iter().zip(&*nekrk).map(|(&a,&b)|f64::min(99.,-10.*f64::log10(error(a,b)))).enumerate().map(|(i,r)| format!("{i}:{r:.0}")).format(" "));
				println!("{} {}", cantera[max.0], nekrk[max.0]);
			}
			max
			//rates.iter().zip(&*cantera).zip(&*nekrk).map(|((&a,&b), &c)| f64::max(num::relative_error(a,b), num::relative_error(b,c))).reduce(f64::max).unwrap()*/
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
		assert!(e.1 < 1e-4, "T {temperature} P {pressure} E {e:?} {max} {}", ""/*species_names[e.0]*/);
		let e = e.1;
		if e > max { max = e; println!("{i} {e:.0e}"); }
	}
	Ok(())
}
