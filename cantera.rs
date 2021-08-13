#![allow(non_snake_case,non_upper_case_globals,mixed_script_confusables)]
#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds,unboxed_closures,fn_traits,trait_alias,iter_zip)]
mod yaml; mod device;
use {anyhow::Result, std::iter::zip, num::sq, iter::{Dot, map}, itertools::Itertools, device::*};
fn main() -> Result<()> {
	let path = std::env::args().skip(1).next().unwrap_or("gri30".to_string());
	let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/share/cantera/data/{path}.yaml") };
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path).expect(&path))?)?;
	let model = yaml::parse(&model);
	use combustion::*;
	let (ref species_names, ref species@Species{ref molar_mass, ref thermodynamics, ..}, _active, _reactions, ref state@State{ref amounts, temperature, pressure_R, ..}) = new(&model);

	use std::os::raw::c_char;
	#[link(name = "cantera")] extern "C" {
		fn thermo_newFromFile(file_name: *const c_char, phase_name: *const c_char) -> i32;
		fn thermo_nSpecies(n: i32) -> usize;
		fn thermo_setTemperature(n: i32, t: f64) -> i32;
		fn thermo_setMoleFractions(n: i32, len: usize, x: *const f64, norm: i32) -> i32;
		fn thermo_getSpeciesName(n: i32, m: usize, len: usize, buffer: *mut c_char) -> i32;
		fn thermo_setPressure(n: i32, p: f64) -> i32;
		fn kin_newFromFile(file_name: *const c_char, phase_name: *const c_char, reactingPhase: i32, neighbor0: i32, neighbor1: i32, neighbor2: i32, neighbor3: i32) -> i32;
	}
	let file = std::ffi::CString::new(path.as_str())?;
	let phase_name_cstr_ptr = std::ffi::CStr::from_bytes_with_nul(if path.contains("gri30.") { b"gri30\0" } else { b"gas\0" })?.as_ptr();
	let phase = unsafe{thermo_newFromFile(file.as_c_str().as_ptr(), phase_name_cstr_ptr)};
	let cantera_species_names = map(0..species.len(), |k| {
		let mut specie = [0; 8];
		unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
		unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
	});
	assert!(unsafe{thermo_nSpecies(phase)} == species.len());
	let _kinetics = unsafe{kin_newFromFile(file.as_c_str().as_ptr(), phase_name_cstr_ptr, phase, 0, 0, 0, 0)};

	#[cfg(not(feature="f32"))] type T = f64;
	#[cfg(feature="f32")] type T = f32;
	#[cfg(not(feature="transport"))] let rates = {
		let rates = reaction::rates(&species.thermodynamics, &reactions);
		with_repetitive_input(assemble::<T>(rates, 1), 1)
	};

	#[cfg(feature="transport")] let transport = {
		eprintln!("Fit");
		let total_amount = amounts.iter().sum::<f64>();
		let mole_fractions = map(&**amounts, |n| n/total_amount);
		let length = 1.;
		let velocity = 1.;
		let time = length / velocity;
		let concentration = pressure_R / temperature;
		let mean_molar_mass = zip(&**molar_mass, &*mole_fractions).map(|(m,x)| m*x).sum::<f64>();
		let density = concentration * mean_molar_mass;
		let Vviscosity = f64::sqrt(density * time) * velocity;
		let mean_molar_heat_capacity_at_CP_R:f64 = thermodynamics.iter().map(|a| a.molar_heat_capacity_at_constant_pressure_R(temperature)).dot(mole_fractions);
		let R = kB*NA;
		let thermal_conductivity0 = mean_molar_heat_capacity_at_CP_R * R / mean_molar_mass * density * length * velocity;
		let mixture_diffusion = sq(length) / time;
		let transport = transport::properties::<5>(&species, temperature, Vviscosity, thermal_conductivity0, density*mixture_diffusion);
		let transport = with_repetitive_input(assemble::<T>(transport, 1), 1);
		{let temperature0 = temperature;
		move |total_amount: T, temperature: T, nonbulk_amounts: &[T]| -> (T, T, Box<[T]>) {
			let_!{ [thermal_conductivity, viscosity, density_mixture_diffusion @ ..] = &*transport(&[], &([&[total_amount, temperature/temperature0], &*nonbulk_amounts].concat())).unwrap() => {
				(thermal_conductivity0*thermal_conductivity, sq(Vviscosity)*viscosity, map(density_mixture_diffusion, |D| density*mixture_diffusion*D))
			}}
		}}
	};

	let ε= 2e-6;
	#[cfg(feature="f32")] let ε = f64::max(ε, 6e-4);
	#[cfg(feature="ir")] let ε = f64::max(ε, 5e-2);

	let test = move |state: &State| -> Result<_> {
		let State{temperature, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();

		let cantera_order = |o: &[f64]| (0..o.len()).map(|i| o[species_names.iter().position(|&s| s==cantera_species_names[i]).unwrap()]).collect::<Box<_>>();
		let amount_fractions = map(&**amounts, |n| n/total_amount);
		unsafe{thermo_setMoleFractions(phase, amount_fractions.len(), cantera_order(&amount_fractions).as_ptr(), 1)}; // /!\ Needs to be set before pressure
		unsafe{thermo_setTemperature(phase, *temperature)};
		unsafe{thermo_setPressure(phase, pressure_R * (kB*NA))}; // /!\ Needs to be set after mole fractions
		let ref nonbulk_amounts = map(&amounts[0..amounts.len()-1], |&n| n as _);
		#[cfg(feature="transport")] {
			#[link(name = "cantera")] extern "C" {
				fn trans_newDefault(th: i32, loglevel: i32) -> i32;
				fn trans_thermalConductivity(n: i32) -> f64;
				fn trans_viscosity(n: i32) -> f64;
				fn trans_getMixDiffCoeffs(n: i32, ldt: i32, dt: *mut f64) -> i32;
			}
			let cantera_transport = unsafe{trans_newDefault(phase, 0)};
			let cantera_thermal_conductivity  = unsafe{trans_thermalConductivity(cantera_transport)};
			let cantera_viscosity = unsafe{trans_viscosity(cantera_transport)};
			let ref cantera_mixture_diffusion = {
				let mut array = vec![0.; species.len()];
				assert!(unsafe{trans_getMixDiffCoeffs(cantera_transport, array.len() as i32, array.as_mut_ptr())} == 0);
				let order = |o: &[f64]| map(&**species_names, |specie| o[cantera_species_names.iter().position(|s| s==specie).unwrap()]);
				order(&map(array, |d| d/*/(pressure_R*NA*kB)*/))
			};
			eprintln!("Cantera");
			eprintln!("λ: {cantera_thermal_conductivity:.4}, μ: {cantera_viscosity:.4e}, D: {:.4e}", cantera_mixture_diffusion.iter().format(" "));
			let State{temperature, amounts, ..} = state;
			let temperature = *temperature as _;
			let total_amount = amounts.iter().sum::<f64>() as _;
			eprintln!("Evaluate");
			let (thermal_conductivity, viscosity, density_mixture_diffusion) = transport(total_amount, temperature, nonbulk_amounts);
			let concentration = pressure_R / temperature;
			let mole_fractions = map(&**amounts, |n| n/total_amount);
			let mean_molar_mass = zip(&**molar_mass, &*mole_fractions).map(|(m,x)| m*x).sum::<f64>();
			let density = concentration * mean_molar_mass;
			let mixture_diffusion = map(&*density_mixture_diffusion, |ρD| ρD/density);
			eprintln!("λ: {thermal_conductivity:.4}, μ: {viscosity:.4e}, D: {:.4e}", mixture_diffusion.iter().format(" "));
			print!("λ: {:.0e}, ", num::relative_error(thermal_conductivity as _, cantera_thermal_conductivity));
			print!("μ: {:.0e}, ", num::relative_error(viscosity as _, cantera_viscosity));
			let e = mixture_diffusion.iter().zip(cantera_mixture_diffusion.iter()).map(|(&a,&b)| num::relative_error(a as _, b)).reduce(f64::max).unwrap();
			println!("D: {e:.0e}");
			Ok((0, 0.))
		}
		#[cfg(not(feature="transport"))] {
			let mut cantera_species_rates = vec![0.; species.len()];
			#[link(name = "cantera")] extern "C" { fn kin_getNetProductionRates(n: i32, len: usize, w_dot: *mut f64) -> i32; }
			unsafe{kin_getNetProductionRates(kinetics, cantera_species_rates.len(), cantera_species_rates.as_mut_ptr())};
			let active = _active;
			let order = |o: &[f64]| map(&species_names[0..active], |specie| o[cantera_species_names.iter().position(|s| s==specie).unwrap()]);
			let cantera_species_rates = order(&map(cantera_species_rates, |c| c*1000.)); // kmol -> mol
			assert!(cantera_species_rates.len() == active);
			// Rates
			let error = |a,b| {
				let m=f64::abs(f64::max(a,b));
				if m < 4e3 { return 0.; }
				let threshold=1e6; (if m < threshold { m/threshold } else { 1. }) * num::relative_error(a,b)
			};
			//let error = |a,b| num::relative_error(a,b);
			let_!{ [_energy_rate_RT, rates @ ..] = &*rates(&[*pressure_R as _], &([&[total_amount as _, *temperature as _], &**nonbulk_amounts].concat()))? => {
			assert!(rates.len() == active, "{}", rates.len());
			if true {
				let (k, e)= rates.iter().zip(&*cantera_rates).map(|(&a,&b)| error(a as _,b)).enumerate().reduce(|a,b| if a.1 > b.1 { a } else { b }).unwrap();
				if e > ε {
					/*println!("{:>10}", species_names.iter().format(" "));
					println!("{:10.0}", cantera_rates.iter().format(" "));
					println!("{:10.0}", rates.iter().format(" "));*/
					println!("{}: {:.0} {:.0}", species_names[k], rates[k], cantera_rates[k]);
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
				use iter::list;
				let lines = list(stdout.split("\n"));
				let ref nekrk_species_names = list(lines[0].trim().split(" "));
				assert!(species_names == nekrk_species_names, "\n{species_names:?}\n{nekrk_species_names:?}");
				let nekrk = map(lines[1].trim().split(" "), |r| r.trim().parse().unwrap());
				assert!(nekrk.len() == nekrk_species_names.len()-1, "{} {}", nekrk.len(), nekrk_species_names.len());
				let (k, e) = cantera_rates.iter().zip(&*nekrk).map(|(&a,&b)| error(a,b)).enumerate().reduce(|a,b| if a.1 > b.1 { a } else { b }).unwrap();
				if e > ε {
					/*println!("{:.0}", cantera_rates.iter().format(" "));
					println!("{:.0}", nekrk.iter().format(" "));
					println!("{:>6}", cantera_rates.iter().zip(&*nekrk).map(|(&a,&b)|f64::min(99.,-10.*f64::log10(error(a,b)))).enumerate().map(|(i,r)| format!("{i}:{r:.0}")).format(" "));*/
					println!("{} {}", cantera_rates[k], nekrk[k]);
				}
				(k, e)
			}}}
			Ok((k, e))
		}
	};
	if false {
		let mut random= rand::thread_rng();
		let mut max = 0.;
		for i in 0.. {
			assert!(state.volume == 1.);
			let volume = state.volume;
			use rand::Rng;
			let temperature = random.gen_range(1000. .. 2800.);
			let pressure = random.gen_range(0. .. 1e5);
			let pressure_R = pressure/(kB*NA); //random.gen_range(0. .. 10e6)/(kB*NA);
			let total_amount = pressure_R * volume / temperature;
			let amount_proportions = map(0..species.len(), |_| random.gen());
			let amounts = map(&*amount_proportions, |amount_proportion| total_amount * amount_proportion/amount_proportions.iter().sum::<f64>());
			let amount_fractions = map(&*amounts, |n| n/total_amount);
			let (k, e) = test(&State{temperature, pressure_R, volume, amounts})?;
			let k = species_names[k];
			assert!(e < ε, "{k}: ε:{e:.0e} T:{temperature:.0}K P:{pressure:.0}Pa X:[{amount_fractions}] (was:{max:.0e})",
				amount_fractions=amount_fractions.iter().format_with(", ",|e,f| f(&format_args!("{:.0}%", e*100.))));
			if e > max { max = e; println!("{i} {e:.0e}"); } else if i%1000==0 { println!("{i}") }
		}
	} else {
			let (_, e) = test(&state)?;
			assert!(e < 1e-7, "{e:e}");
	}
	Ok(())
}
