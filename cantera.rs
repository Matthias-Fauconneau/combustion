#![allow(non_snake_case,non_upper_case_globals,mixed_script_confusables)]
#![feature(format_args_capture,default_free_fn,associated_type_bounds,unboxed_closures,fn_traits,trait_alias,iter_zip,let_else)]
mod yaml; mod device;
use {std::iter::zip, anyhow::Result, iter::{map, Dot}, itertools::Itertools, combustion::*, device::*};
fn main() -> Result<()> {
	let path = std::env::args().skip(1).next().unwrap_or("LiDryer".to_string());
	let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/local/share/cantera/data/{path}.yaml") };
	let yaml = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path).expect(&path))?)?;
	let model = yaml::parse(&yaml);
	let (ref species_names, ref species@Species{ref molar_mass, ref thermodynamics, ..}, _active, _reactions, /*ref state@State{ref amounts, temperature, pressure_R, ..}*/..) = new(&model);
	let state = State{
		pressure_R: 101325.0 / R,
		volume: 1.,
		temperature: 1556.971923828125,
		amounts: vec![0.09260429441928864, 0.052496496587991717, 0.7380790114402771, 0.07698292285203934, 0.012342223897576332, 0.02748161181807518, 0.000013450758160615806, 6.392859575043985e-10, 1.880952380952381].into_boxed_slice()
	};
	let ref state@State{ref amounts, temperature, pressure_R, ..} = state;
	let _ = (amounts, temperature, pressure_R);
	eprintln!("{}", zip(&**amounts,&**species_names).filter(|(&x,_)| x != 0.).map(|(_,name)| name).format(" "));

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
		let mut specie = [0; 32];
		unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
		unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
	});
	assert!(unsafe{thermo_nSpecies(phase)} == species.len());
	let _kinetics = unsafe{kin_newFromFile(file.as_c_str().as_ptr(), phase_name_cstr_ptr, phase, 0, 0, 0, 0)};

	#[cfg(not(feature="transport"))] let rates = {
		let rates = reaction::rates(molar_mass, thermodynamics, &reactions, species_names);
		with_repetitive_input(assemble::<f32>(rates, 1), 1)
	};

	#[cfg(feature="transport")] let transport = {
		//std::iter::zip, iter::{Dot,
		eprintln!("Fit");
		let (temperature0, viscosity0, thermal_conductivity0) = if true {
			let total_amount = amounts.iter().sum::<f64>();
			let mole_fractions = map(&**amounts, |n| n/total_amount);
			let diffusivity = 1.;
			let concentration = pressure_R / temperature;
			let mean_molar_mass = zip(&**molar_mass, &*mole_fractions).map(|(m,x)| m*x).sum::<f64>();
			let density = concentration * mean_molar_mass;
			let viscosity = density * diffusivity;
			let mean_molar_heat_capacity_at_CP_R:f64 = thermodynamics.iter().map(|a| a.molar_heat_capacity_at_constant_pressure_R(temperature)).dot(mole_fractions);
			let thermal_conductivity = mean_molar_heat_capacity_at_CP_R * R / mean_molar_mass * viscosity;
			(temperature, viscosity, thermal_conductivity)
		} else { (1., 1., 1.) };
		let transport = transport::properties::<5>(&species, temperature0, viscosity0, thermal_conductivity0);
		let transport = with_repetitive_input(assemble(transport, 1), 1);
		move |total_amount:f64, temperature:f64, nonbulk_amounts:&[f64]| -> (f64, f64, Box<[f64]>) {
			let nonbulk_amounts = map(&*nonbulk_amounts, |&n| n as f32);
			let [thermal_conductivity, viscosity, density_diffusion @ ..] = &*transport(&[], &([&[total_amount as f32, (temperature/temperature0) as f32], &*nonbulk_amounts].concat())).unwrap() else { panic!() };
			(thermal_conductivity0*(*thermal_conductivity as f64), viscosity0*(*viscosity as f64), map(density_diffusion, |&D:&f32| viscosity0*(D as f64)))
		}
	};

	let ε= 1e-3;
	#[cfg(feature="ir")] let ε = f64::max(ε, 1e-2);

	let test = move |state: &State| -> Result<_> {
		let State{temperature, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();

		let cantera_order = |o: &[f64]| (0..o.len()).map(|i| o[species_names.iter().position(|&s| s==cantera_species_names[i]).expect(&cantera_species_names[i])]).collect::<Box<_>>();
		let amount_fractions = map(&**amounts, |n| n/total_amount);
		unsafe{thermo_setMoleFractions(phase, amount_fractions.len(), cantera_order(&amount_fractions).as_ptr(), 1)}; // /!\ Needs to be set before pressure
		unsafe{thermo_setTemperature(phase, *temperature)};
		unsafe{thermo_setPressure(phase, pressure_R * R)}; // /!\ Needs to be set after mole fractions

		let nonbulk_amounts = &amounts[0..amounts.len()-1];
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
			let ref cantera_diffusion = {
				let mut array = vec![0.; species.len()];
				assert!(unsafe{trans_getMixDiffCoeffs(cantera_transport, array.len() as i32, array.as_mut_ptr())} == 0);
				let order = |o: &[f64]| map(&**species_names, |specie| o[cantera_species_names.iter().position(|s| s==specie).unwrap()]);
				order(&map(array, |d| d/*/(pressure_R*NA*kB)*/))
			};
			eprintln!("Cantera");
			eprintln!("λ: {cantera_thermal_conductivity:.4}, μ: {cantera_viscosity:.4e}, D: {:.4e}", cantera_diffusion.iter().format(" "));
			let State{temperature, amounts, ..} = state;
			let temperature = *temperature;
			let total_amount = amounts.iter().sum::<f64>();
			eprintln!("Evaluate");
			let (thermal_conductivity, viscosity, density_diffusion) = transport(total_amount, temperature, nonbulk_amounts);
			let concentration = pressure_R / temperature;
			let mole_fractions = map(&**amounts, |n| n/total_amount);
			let mean_molar_mass = zip(&**molar_mass, &*mole_fractions).map(|(m,x)| m*x).sum::<f64>();
			let density = concentration * mean_molar_mass;
			let diffusion = map(&*density_diffusion, |ρD| ρD/density);
			eprintln!("λ: {thermal_conductivity:.4}, μ: {viscosity:.4e}, D: {:.4e}", diffusion.iter().format(" "));
			print!("λ: {:.0e}, ", num::relative_error(thermal_conductivity as _, cantera_thermal_conductivity));
			print!("μ: {:.0e}, ", num::relative_error(viscosity as _, cantera_viscosity));
			let e = diffusion.iter().zip(cantera_diffusion.iter()).map(|(&a,&b)| num::relative_error(a as _, b)).reduce(f64::max).unwrap();
			println!("D: {e:.0e}");
			Ok((0, 0.))
		}
		#[cfg(not(feature="transport"))] {
			let mut cantera_species_rates = vec![0.; species.len()];
			#[link(name = "cantera")] extern "C" { fn kin_getNetProductionRates(n: i32, len: usize, w_dot: *mut f64) -> i32; }
			unsafe{kin_getNetProductionRates(kinetics, cantera_species_rates.len(), cantera_species_rates.as_mut_ptr())};
			let active = _active;
			let order = |o: &[f64]| map(&species_names[0..active], |specie| o[cantera_species_names.iter().position(|s| s==specie).expect(specie)]);
			let cantera_species_rates = order(&map(cantera_species_rates, |c| c*1000.)); // kmol -> mol
			assert!(cantera_species_rates.len() == active);
			// Rates
			let error = |x,r| f64::abs(x-r) / f64::max(r, 1e6) /*{
				let m=f64::abs(f64::max(a,b));
				if m < 4e3 { return 0.; }
				let threshold=1e6; (if m < threshold { m/threshold } else { 1. }) * num::relative_error(a,b)
			}*/;
			let [_dtT, _dtV, rates @ ..] = &*rates(&[*pressure_R as _, (1./pressure_R) as _], &([&[total_amount as f32, *temperature as _] as &[_], &*nonbulk_amounts].concat()))? else { panic!() };
			assert!(rates.len() == active, "{}", rates.len());
			if true {
				let (k, e)= rates.iter().zip(&*cantera_species_rates).map(|(&x,&r)| error(x as _,r)).enumerate().reduce(|a,b| if a.1 > b.1 { a } else { b }).unwrap();
				if e > ε {
					/*println!("{:>10}", species_names.iter().format(" "));
					println!("{:10.0}", cantera_species_rates.iter().format(" "));
					println!("{:10.0}", rates.iter().format(" "));*/
					println!("{}: {:.0} {:.0}", species_names[k], rates[k], cantera_species_rates[k]);
				}
				Ok((k, e))
			}
			else {
				//assert!(std::process::Command::new("make").current_dir("../nekRK/build").arg("-j").status()?.success());
				let code = std::path::Path::new("/var/tmp/main.c");
				if !code.exists() || std::fs::metadata("../nekRK/main.py")?.modified()? > code.metadata()?.modified()? {
					std::fs::write("/var/tmp/main.c", std::process::Command::new("../nekRK/main.py").arg(&path).output()?.stdout)?;
				}
				let nekrk = std::process::Command::new("../nekRK/build/main").current_dir("../nekRK")
					.args([code.to_str().unwrap(),"Serial","1","1","0",&(pressure_R*R).to_string(),&temperature.to_string(),&amount_fractions.iter().format(" ").to_string()]).output()?;
				assert!(nekrk.stderr.is_empty(), "{}", std::str::from_utf8(&nekrk.stderr)?);
				let stdout = nekrk.stdout;
				let stdout = std::str::from_utf8(&stdout)?;
				use iter::list;
				let lines = list(stdout.split("\n"));
				let ref nekrk_species_names = list(lines[0].trim().split(" "));
				assert!(species_names == nekrk_species_names, "\n{species_names:?}\n{nekrk_species_names:?}");
				let nekrk = map(lines[1].trim().split(" "), |r| r.trim().parse().unwrap());
				assert!(nekrk.len() == nekrk_species_names.len()-1, "{} {}", nekrk.len(), nekrk_species_names.len());
				let (k, e) = cantera_species_rates.iter().zip(&*nekrk).map(|(&a,&b)| error(a,b)).enumerate().reduce(|a,b| if a.1 > b.1 { a } else { b }).unwrap();
				if e > ε {
					/*println!("{:.0}", cantera_species_rates.iter().format(" "));
					println!("{:.0}", nekrk.iter().format(" "));
					println!("{:>6}", cantera_species_rates.iter().zip(&*nekrk).map(|(&a,&b)|f64::min(99.,-10.*f64::log10(error(a,b)))).enumerate().map(|(i,r)| format!("{i}:{r:.0}")).format(" "));*/
					println!("{} {}", cantera_species_rates[k], nekrk[k]);
				}
				Ok((k, e))
			}
		}
	};
	if false {
		let mut random= rand::thread_rng();
		let mut max = 0.;
		let start = std::time::Instant::now();
		for i in 0.. {
			assert!(state.volume == 1.);
			let volume = state.volume;
			use rand::Rng;
			let temperature = random.gen_range(1000. .. 2700.);
			let pressure = random.gen_range(1000. .. 1e5);
			let pressure_R = pressure/R; //random.gen_range(0. .. 10e6)/R;
			let total_amount = pressure_R * volume / temperature;
			let amount_proportions = map(&**amounts, |&x| random.gen::<f64>() * if x != 0. { 1. } else { 1e-17 });
			let amounts = map(&*amount_proportions, |amount_proportion| total_amount * amount_proportion/amount_proportions.iter().sum::<f64>());
			let amount_fractions = map(&*amounts, |n| n/total_amount);
			//println!("T:{temperature:.0}K P:{pressure:.0}Pa X:[{amount_fractions}]", amount_fractions=amount_fractions.iter().format_with(", ",|e,f| f(&format_args!("{:.0}%", e*100.))));
			let (k, e) = test(&State{temperature, pressure_R, volume, amounts})?;
			let k = species_names[k];
			assert!(e < ε, "{k}: ε:{e:.0e} T:{temperature:.0}K P:{pressure:.0}Pa X:[{amount_fractions}] (was:{max:.0e})",
				amount_fractions=zip(&**species_names, amount_fractions.iter()).format_with(", ",|(specie,e),f| f(&format_args!("{specie}:{:.0}%", e*100.))));
			if e > max { max = e; println!("{i} {e:.0e}"); } else if i%100000==0 { println!("{i} {max:.0e} {:.0}K/s", (i as f64/1000.)/start.elapsed().as_secs_f64()) }
		}
	} else {
		println!("{state:?}");
		let (_, e) = test(&state)?;
		println!("{e:e}");
		assert!(e < 2e-3, "{e:e}");
	}
	Ok(())
}
