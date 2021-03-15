pub fn promote(v: &[f32]) -> Box<[f64]> { v.iter().map(|&v| v as f64).collect() }
use super::*;
use itertools::Itertools;

#[throws] pub fn check(model: Model, Simulation{species_names, /*time_step,*/ state, ..}: &Simulation) {
	let len = model.len();
	let file = std::ffi::CStr::from_bytes_with_nul(b"gri30.yaml\0").unwrap().as_ptr();
	let name = std::ffi::CStr::from_bytes_with_nul(b"gri30\0").unwrap().as_ptr();
	let phase = unsafe{thermo_newFromFile(file, name)};
	assert!(unsafe{thermo_nSpecies(phase)} == len);
	{
		let cantera_species_name = (0..len).map(|k| {
			let mut specie = [0; 8];
			unsafe{thermo_getSpeciesName(phase, k, specie.len(), specie.as_mut_ptr())};
			unsafe{std::ffi::CStr::from_ptr(specie.as_ptr()).to_str().unwrap().to_owned()}
		}).collect::<Box<_>>();
		assert_eq!(&cantera_species_name.iter().map(String::as_str).collect::<Box<_>>(), species_names);
	}
	//let time = 0.;
	const CONSTANT : Property = {use Property::*; Volume};
	let volume = state.volume;
	//let mut cvode = cvode::CVODE::new(promote(((&state).into():StateVector<CONSTANT>).0.deref()).deref());
	loop {
		//let next_time = time + *time_step;
		let ref cantera_rate = {
			assert!(state.amounts.len() == len);
			unsafe{thermo_setMoleFractions(phase, state.amounts.len(), state.amounts.as_ptr(), 1)}; // /!\ Needs to be set before pressure
			unsafe{thermo_setTemperature(phase, state.temperature)};
			unsafe{thermo_setPressure(phase, state.pressure * NA)}; // /!\ Needs to be set after mole fractions
			let kinetics = unsafe{kin_newFromFile(file, name, phase, 0, 0, 0, 0)};
			let mut net_productions_rates = vec![0.; len];
			unsafe{kin_getNetProductionRates(kinetics, len, net_productions_rates.as_mut_ptr())};
			net_productions_rates.iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>() // kmol/m^3/s => mol/s [1m^3]
		};
		let (_, _, rate) = model.rate();
		let mut derivative = /*Derivative*/StateVector::<CONSTANT>(std::iter::repeat(0.).take(len).collect());
		let rate = {rate(state.constant(), &state.into(), &mut derivative); derivative.0};

		fn table(labels: &[&str], a: &[f64], b: &[f64]) -> Box<[([String; 4], usize)]> {
			labels.iter().zip(a.iter().zip(b)).filter(|(_,(&a,&b))| a.abs() > 1e-29 || b.abs() > 1e-29).map(|(&header,(&a,&b))| {
				fn to_string(v: f64) -> String { if v == 0. { "0".to_owned() } else { format!("{:.0e}", v) } }
				let column = [header.to_owned(), to_string(a), to_string(b), to_string(num::relative_error(a,b))];
				let width = column.iter().map(|s| s.len()).max().unwrap();
				(column, width)
			}).collect()
		}
		fn print<const R: usize>(table: &[([String; R], usize)]) {
			for row in 0..R { println!("{}", table.iter().format_with(" ", |(c,width), f| f(&format_args!("{:width$}", c[row], width=width)))); }
		}

		let rate = &rate[2..];
		let rate = promote(&rate);
		let ref rate = rate.iter().map(|dtn| dtn/volume).collect::<Box<_>>();
		let ref rate = rate.iter().map(|&dtn| if dtn.abs() < 1e-29 { 0. } else { dtn }).collect::<Box<_>>(); // Correct max relative error
		let len = species_names.len();
		print(&table(&species_names[..len-1], rate, &cantera_rate[..len-1]));
		fn absolute_error(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(&a,&b)| f64::abs(a-b)).reduce(f64::max).unwrap() }
		fn relative_error(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(&a,&b)| num::relative_error(a,b)).reduce(f64::max).unwrap() }
		{
			let abs = absolute_error(rate, cantera_rate);
			let rel = relative_error(rate, cantera_rate);
			println!("rate {:e} {:e}", abs, rel);
			assert!(abs < 1e-8 && rel < 4e-4, "rate {:e} {:e}", abs, rel);
		}

		/*while time < next_time {
			let (next_time, next_state) = cvode.step(move |u| system.rate_and_jacobian::<CONSTANT>(state.constant::<CONSTANT>(), &State(*u)).map(|(rate, /*jacobian*/)| rate.0), next_time, &((&state).into():reaction::State<CONSTANT,S>)); //dbg!(time);
			(time, state) = (next_time, State::new(state.amounts.iter().sum(), state.constant::<CONSTANT>(), &reaction::State::<CONSTANT, S>(next_state)))
		}
		assert_eq!(time, next_time);
		println!("t {}", time);*/

		/*{
			println!("T {} {} {:e}", state.temperature, cantera_state.temperature, num::relative_error(state.temperature, cantera_state.temperature));
			print(&table(species_names, &state.amounts, &cantera_state.amounts));
			{
				let abs = absolute_error(&state.amounts, &cantera_state.amounts);
				let rel = relative_error(&state.amounts, &cantera_state.amounts);
				println!("state {:e} {:e}", abs, rel);
				assert!(abs < 1e-8 || rel < 0., "state {:e} {:e}", abs, rel);
			}
		}*/

		//state = *cantera_state; // Check rates along cantera trajectory
		break;
	}
}
