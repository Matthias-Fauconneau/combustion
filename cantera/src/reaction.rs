pub fn map<T, U, F: Fn(&T)->U>(v: &[T], f: F) -> Box<[U]> { v.iter().map(f).collect() }
pub fn promote(v: &[f32]) -> Box<[f64]> { map(v, |&v| v as f64) }
pub fn demote(v: &[f64]) -> Box<[f32]> { map(v, |&v| v as f32) }
fn implicit(u: &[f64]) -> Box<[f64]> { [u[0]].iter().chain(&u[2..]).copied().collect() } // one of P or V imply the other using ideal gas law
fn explicit(total_amount: f64, pressure: f64, u: &[f64]) -> Box<[f64]> { // Reconstructs P|V using ideal gas law
	let temperature = u[0];
	let volume = K * total_amount / pressure * temperature;
	[temperature, volume].iter().chain(&u[1..]).copied().collect()
}

use super::{*, Property::*};

#[throws] pub fn check(model: Model, Simulation{species_names, time_step, state, ..}: &Simulation) {
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
	let total_amount = state.amounts.iter().sum();
	let volume = state.volume;
	let constant = state.constant();
	let mut state = implicit(&promote(&(state.into():StateVector<{Volume}>)));
	let mut cvode = cvode::CVODE::new(&state);
	let mut time = 0.;
	let derivative = /*Derivative*/StateVector(std::iter::repeat(0.).take(2+len-1).collect());
	let derivative = std::cell::Cell::new(derivative);
	fn get_mut<T: Default, U>(cell: &std::cell::Cell<T>, f: impl FnOnce(&mut T) -> U) -> U {
		let mut value = cell.take();
		let result = f(&mut value);
		cell.set(value);
		result
	}
	let (_, rate) = model.rate();
	/*let ref derivative = move |u| get_mut(derivative, |derivative| {
		let pressure = constant.0 as f64;
		rate(constant, &StateVector(map(&explicit(total_amount, pressure, u), |&v| (v as f32).max(0.))), &mut derivative);
		Some(implicit(&promote(&derivative.0)))
	});*/
	// CVODE shim
	struct Derivative<Rate: crate::Rate<CONSTANT>, const CONSTANT: Property> {
		rate: Rate,
		constant: Constant<CONSTANT>,
		total_amount: f64,
		derivative: std::cell::Cell<StateVector::<CONSTANT>>
	}
	impl<Rate: crate::Rate<CONSTANT>, const CONSTANT: Property> FnOnce<(&[f64],)> for Derivative<Rate, CONSTANT> {
		type Output = Option<Box<[f64]>>;
		extern "rust-call" fn call_once(mut self, args: (&[f64],)) -> Self::Output { self.call_mut(args) }
	}
	impl<Rate: crate::Rate<CONSTANT>, const CONSTANT: Property> FnMut<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call_mut(&mut self, args: (&[f64],)) -> Self::Output { self.call(args) }
	}
	impl<Rate: crate::Rate<CONSTANT>, const CONSTANT: Property> Fn<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call(&self, (u,): (&[f64],)) -> Self::Output {
			let Self{rate, constant, total_amount, derivative} = self;
			get_mut(derivative, |mut derivative| {
				let pressure = constant.0 as f64;
				rate(*constant, &StateVector(map(&explicit(*total_amount, pressure, u), |&v| (v as f32).max(0.))), &mut derivative);
				Some(implicit(&promote(&derivative.0)))
			})
		}
	}
	let ref derivative = Derivative{
		rate,
		constant,
		total_amount,
		derivative
	};

	let mut last_time = time;
	while std::hint::black_box(true) {
		let ref state_vector = StateVector(demote(&explicit(total_amount, constant.0 as f64, &state)));
		let ref cantera_rate = {
			let state = State::new(total_amount, constant, state_vector);
			assert!(state.amounts.len() == len);
			unsafe{thermo_setMoleFractions(phase, state.amounts.len(), state.amounts.as_ptr(), 1)}; // /!\ Needs to be set before pressure
			unsafe{thermo_setTemperature(phase, state.temperature)};
			unsafe{thermo_setPressure(phase, state.pressure * NA)}; // /!\ Needs to be set after mole fractions
			let kinetics = unsafe{kin_newFromFile(file, name, phase, 0, 0, 0, 0)};
			let mut net_productions_rates = vec![0.; len];
			unsafe{kin_getNetProductionRates(kinetics, len, net_productions_rates.as_mut_ptr())};
			net_productions_rates.iter().map(|c| c*1000.).take(len-1).collect::<Box<_>>() // kmol -> mol
		};
		let rate = {
			let ref rate = derivative.rate;
			let mut derivative = /*Derivative*/StateVector::<{Volume}>(std::iter::repeat(0.).take(2+len-1).collect());
			rate(constant, state_vector, &mut derivative);
			derivative.0
		};

		let rate = &rate[2..];
		let rate = promote(&rate);
		let ref rate = rate.iter().map(|dtn| dtn/volume).collect::<Box<_>>();
		//let ref rate = rate.iter().map(|&dtn| if dtn.abs() < 1e-29 { 0. } else { dtn }).collect::<Box<_>>(); // Correct max relative error
		fn absolute_error(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(&a,&b)| f64::abs(a-b)).reduce(f64::max).unwrap() }
		fn relative_error(a: &[f64], b: &[f64]) -> f64 {
			a.iter().zip(b).map(|(&a,&b)|
				if a.abs() < 1e-3 && b.abs() < 1e-3 { 0. } else { num::relative_error(a,b) }
			).reduce(f64::max).unwrap()
		}
		let abs = absolute_error(rate, cantera_rate);
		let rel = relative_error(rate, cantera_rate);
		//println!("{} {:e} {:e} {:e}", time*1e3, time-last_time, abs, rel);
		last_time = time;
		if !(abs < 1e-4 && rel < 1e-2) {
			fn table(labels: &[&str], a: &[f64], b: &[f64]) -> Box<[([String; 5], usize)]> {
				labels.iter().zip(a.iter().zip(b)).filter(|(_,(&a,&b))| a.abs() > 1e-29 || b.abs() > 1e-29).map(|(&header,(&a,&b))| {
					fn to_string(v: f64) -> String { if v == 0. { "0".to_owned() } else { format!("{:.0e}", v) } }
					let column = [header.to_owned(), to_string(a), to_string(b), to_string(num::abs(a-b)), to_string(num::relative_error(a,b))];
					let width = column.iter().map(|s| s.len()).max().unwrap();
					(column, width)
				}).collect()
			}
			fn print<const R: usize>(table: &[([String; R], usize)]) {
				use itertools::Itertools;
				for row in 0..R { println!("{}", table.iter().format_with(" ", |(c,width), f| f(&format_args!("{:width$}", c[row], width=width)))); }
			}
			print(&table(&species_names[..len-1], rate, &cantera_rate[..len-1]));
		}
		assert!(abs < 1e-4 && rel < 1e-2, "{:e} {:e}", abs, rel);

		/*while time < next_time {
			let (next_time, next_state) = cvode.step(move |u| system.rate_and_jacobian::<CONSTANT>(state.constant::<CONSTANT>(), &State(*u)).map(|(rate, /*jacobian*/)| rate.0), next_time, &((&state).into():reaction::State<CONSTANT,S>)); //dbg!(time);
			(time, state) = (next_time, State::new(state.amounts.iter().sum(), state.constant::<CONSTANT>(), &reaction::State::<CONSTANT, S>(next_state)))
		}
		assert_eq!(time, next_time);
		println!("t {}", time);*/
		(time, state) = {
			let (time, state) = cvode.step(derivative, time+time_step, &state);
			(time, state.to_vec().into_boxed_slice())
		}

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
		//break;
	}
}
