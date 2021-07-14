#![feature(format_args_capture,array_map,in_band_lifetimes,default_free_fn,associated_type_bounds,unboxed_closures,fn_traits)]#![allow(non_snake_case,non_upper_case_globals)]
mod yaml; mod device;
use {anyhow::Result, iter::{list, map}, itertools::Itertools, device::*};
fn main() -> Result<()> {
	let path = std::env::args().skip(1).next().unwrap();
	let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/share/cantera/data/{path}.yaml") };
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path).expect(&path))?)?;
	let model = yaml::parse(&model);
	use combustion::*;
	let (ref species_names, ref species, active, reactions, state) = new(&model);

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
	fn trans_newDefault(th: i32, loglevel: i32) -> i32;
	fn trans_thermalConductivity(n: i32) -> f64;
	fn trans_viscosity(n: i32) -> f64;
	fn trans_getMixDiffCoeffs(n: i32, ldt: i32, dt: *mut f64) -> i32;
	}
	let file = std::ffi::CString::new(path.as_str())?;
	let phase_name_cstr_ptr = std::ffi::CStr::from_bytes_with_nul(if path.contains("gri30") { b"gri30\0" } else { b"gas\0" })?.as_ptr();
	let phase = unsafe{thermo_newFromFile(file.as_c_str().as_ptr(), phase_name_cstr_ptr)};
	let cantera_species_names = map(0..species.len(), |k| {
		let mut specie = [0; 8];
		unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
		unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
	});
	assert!(unsafe{thermo_nSpecies(phase)} == species.len());
	let kinetics = unsafe{kin_newFromFile(file.as_c_str().as_ptr(), phase_name_cstr_ptr, phase, 0, 0, 0, 0)};

	let test = move |state: &State| -> Result<_> {
		let State{temperature, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();

		let cantera_order = |o: &[f64]| (0..o.len()).map(|i| o[species_names.iter().position(|&s| s==cantera_species_names[i]).unwrap()]).collect::<Box<_>>();
		let amount_fractions = map(&**amounts, |n| n/total_amount);
		unsafe{thermo_setMoleFractions(phase, amount_fractions.len(), cantera_order(&amount_fractions).as_ptr(), 1)}; // /!\ Needs to be set before pressure
		unsafe{thermo_setTemperature(phase, *temperature)};
		unsafe{thermo_setPressure(phase, pressure_R * (kB*NA))}; // /!\ Needs to be set after mole fractions
		let transport = unsafe{trans_newDefault(phase, 0)};
    let cantera_thermal_conductivity  = unsafe{trans_thermalConductivity(transport)};
    let cantera_viscosity = unsafe{trans_viscosity(transport)};
    let ref cantera_mixture_diffusion_coefficients = {
			let mut array = vec![0.; species.len()];
			assert!(unsafe{trans_getMixDiffCoeffs(transport, array.len() as i32, array.as_mut_ptr())} == 0);
			let order = |o: &[f64]| map(&**species_names, |specie| o[cantera_species_names.iter().position(|s| s==specie).unwrap()]);
			order(&map(array, |d| d/*/(pressure_R*NA*kB)*/))
		};
		eprintln!("Cantera");
		eprintln!("λ: {cantera_thermal_conductivity:.4}, μ: {cantera_viscosity:.4e}, D: {:.4e}", cantera_mixture_diffusion_coefficients.iter().format(" "));

		let ref nonbulk_amounts = map(&amounts[0..amounts.len()-1], |&n| n as _);
		#[cfg(feature="transport")] {
		eprintln!("Fit");
		let transport = /*if false {
			let transport = transport::properties::<5>(&species);
			let transport = with_repetitive_input(assemble(&transport), 1);
			|pressure_R: f64, total_amount: f64, temperature: f64, nonbulk_amounts: &[f64]| -> (f64, f64, Box<[f64]>) {
				let_!{ [thermal_conductivity, viscosity, mixture_diffusion_coefficients @ ..] = &*transport(&[pressure_R as _], &([&[total_amount as _, temperature as _], &*nonbulk_amounts].concat())).unwrap() => {
					(*thermal_conductivity, *viscosity, mixture_diffusion_coefficients.into())
				}}
			}
		} else {*/
			transport::properties_rust::<5>(&species)
		//}
		;
		let State{temperature, pressure_R, amounts, ..} = state;
		let temperature = *temperature;
		let pressure_R = *pressure_R;
		let total_amount = amounts.iter().sum::<f64>();
		eprintln!("Evaluate");
		let (thermal_conductivity, viscosity, mixture_diffusion_coefficients) = transport(pressure_R, total_amount, temperature, nonbulk_amounts);
		eprintln!("λ: {thermal_conductivity:.4}, μ: {viscosity:.4e}, D: {:.4e}", mixture_diffusion_coefficients.iter().format(" "));
    print!("λ: {:.0e}, ", num::relative_error(thermal_conductivity as _, cantera_thermal_conductivity));
    print!("μ: {:.0e}, ", num::relative_error(viscosity as _, cantera_viscosity));
    let e = mixture_diffusion_coefficients.iter().zip(cantera_mixture_diffusion_coefficients.iter()).map(|(&a,&b)| num::relative_error(a as _, b)).reduce(f64::max).unwrap();
    println!("D: {e:.0e}");
		}
		let (k, e) = if false {
			let cantera_rates = if true { // Species
				let mut rates = vec![0.; species.len()];
				unsafe{kin_getNetProductionRates(kinetics, rates.len(), rates.as_mut_ptr())};
				let order = |o: &[f64]| map(&species_names[0..active], |specie| o[cantera_species_names.iter().position(|s| s==specie).unwrap()]);
				let species_rates = order(&map(rates, |c| c*1000.)); // kmol -> mol
				assert!(species_rates.len() == active);
				species_rates
			} else { // Reactions
				let mut rates = vec![0.; reactions.len()];
				unsafe{kin_getNetRatesOfProgress(kinetics, rates.len(), rates.as_mut_ptr())};
				map(rates, |c| c*1000.) // kmol -> mol
			};

		// Rates
		let error = |a,b| {
			let m=f64::abs(f64::max(a,b));
			if m < 1e2 { return 0.; }
			let threshold=1e6; (if m < threshold { m/threshold } else { 1. }) * num::relative_error(a,b)
		};
		//let error = |a,b| num::relative_error(a,b);
		let rates = reaction::rates(&species.thermodynamics, &reactions);
		let rates = with_repetitive_input(assemble(&rates), 1);
		let_!{ [_energy_rate_RT, rates @ ..] = &*rates(&[*pressure_R as _], &([&[total_amount as _, *temperature as _], &**nonbulk_amounts].concat()))? => {
		assert!(rates.len() == active, "{}", rates.len());
		if true {
			let (k, e)= rates.iter().zip(&*cantera_rates).map(|(&a,&b)| error(a as _,b)).enumerate().reduce(|a,b| if a.1 > b.1 { a } else { b }).unwrap();
			if e > 1e-4 {
				println!("{:.0}", cantera_rates.iter().format(" "));
				println!("{}", rates.iter().format(" "));
				println!("{} {}", rates[k], cantera_rates[k]);
			}
			(k, e)
		}
		else {
			//assert!(std::process::Command::new("make").current_dir("../nekRK/build").arg("-j").status()?.success());
			let code = std::path::Path::new("/var/tmp/main.c");
			if !code.exists() || std::fs::metadata("../nekRK/main.py")?.modified()? > code.metadata()?.modified()? {
				std::fs::write("/var/tmp/main.c", std::process::Command::new("../nekRK/main.py").arg(&path).output()?.stdout)?;
			}
			let nekrk = std::process::Command::new("../nekRK/build/main").current_dir("../nekRK")
				.args([code.to_str().unwrap(),"Serial","1","1","0",&(pressure_R*(kB*NA)).to_string(),&temperature.to_string(),&amount_fractions.iter().format(" ").to_string()]).output()?;
			assert!(nekrk.stderr.is_empty(), "{}", std::str::from_utf8(&nekrk.stderr)?);
			let stdout = nekrk.stdout;
			let stdout = std::str::from_utf8(&stdout)?;
			let lines = list(stdout.split("\n"));
			let ref nekrk_species_names = list(lines[0].trim().split(" "));
			assert!(species_names == nekrk_species_names, "\n{species_names:?}\n{nekrk_species_names:?}");
			let nekrk = map(lines[1].trim().split(" "), |r| r.trim().parse().unwrap());
			assert!(nekrk.len() == nekrk_species_names.len()-1, "{} {}", nekrk.len(), nekrk_species_names.len());
			let (k, e) = cantera_rates.iter().zip(&*nekrk).map(|(&a,&b)| error(a,b)).enumerate().reduce(|a,b| if a.1 > b.1 { a } else { b }).unwrap();
			if e > 1e-4 {
				println!("{:.0}", cantera_rates.iter().format(" "));
				println!("{:.0}", nekrk.iter().format(" "));
				println!("{:>6}", cantera_rates.iter().zip(&*nekrk).map(|(&a,&b)|f64::min(99.,-10.*f64::log10(error(a,b)))).enumerate().map(|(i,r)| format!("{i}:{r:.0}")).format(" "));
				println!("{} {}", cantera_rates[k], nekrk[k]);
			}
			(k, e)
		}}}
		} else { (0, 0.) };
		Ok((k, e))
	};
	if false {
		let mut random= rand::thread_rng();
		let mut max = 0.;
		for i in 0.. {
			assert!(state.volume == 1.);
			let volume = state.volume;
			use rand::Rng;
			let temperature = random.gen_range(1000. .. 3500.);
			let pressure = random.gen_range(0. .. 1e5);
			let pressure_R = pressure/(kB*NA); //random.gen_range(0. .. 10e6)/(kB*NA);
			let total_amount = pressure_R * volume / temperature;
			let amount_proportions = map(0..species.len(), |_| random.gen());
			let amounts = map(&*amount_proportions, |amount_proportion| total_amount * amount_proportion/amount_proportions.iter().sum::<f64>());
			let amount_fractions = map(&*amounts, |n| n/total_amount);
			let (k, e) = test(&State{temperature, pressure_R, volume, amounts})?;
			assert!(e < 1e-13, "T:{temperature} P:{pressure} X:{amount_fractions:?} E:{}:{:e}(vs{max:e})", species_names[k], e);
			if e > max { max = e; println!("{i} {e:.0e}"); } else if i%100==0 { println!("{i}") }
		}
	} else {
			let (_, e) = test(&state)?;
			assert!(e < 1e-13);
	}
	Ok(())
}
