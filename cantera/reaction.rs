fn get_mut<T: Default, U>(cell: &std::cell::Cell<T>, f: impl FnOnce(&mut T) -> U) -> U {
	let mut value = cell.take();
	let result = f(&mut value);
	cell.set(value);
	result
}

fn implicit(u: &[f64]) -> Box<[f64]> { [u[0]].iter().chain(&u[2..]).copied().collect() } // one of P or V imply the other using ideal gas law
fn explicit(total_amount: f64, pressure_R: f64, u: &[f64]) -> Box<[f64]> { // Reconstructs P|V using ideal gas law
	let temperature = u[0];
	let volume = total_amount / pressure_R * temperature;
	[temperature, volume].iter().chain(&u[1..]).copied().collect()
}

use combustion::{*, reaction::{*, Property::*}};
pub fn check(model: &model::Model, state: &State) {
	let (ref species_names, ref species) = Species::new(&model.species);
	let len = species.len();
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(species_names, r));
	assert!(state.amounts.len() == len && !state.amounts.iter().any(|&n| n<0.));
	let total_amount = state.amounts.iter().sum();
	let constant = state.constant();
	let state : StateVector<{Volume}> = state.into();
	let state = implicit(&state);
	let mut cvode = cvode::CVODE::new(/*relative_tolerance:*/ 1e-4, /*absolute_tolerance:*/ 1e-7, &state);
	let derivative = /*Derivative*/StateVector::<{Pressure}>(std::iter::repeat(0.).take(2+len-1).collect());
	let derivative = std::cell::Cell::new(derivative);
	let rate = rate_function(species, &*reactions);
	// CVODE shim
	struct Derivative<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> {
		rate: Rate,
		constant: Constant<CONSTANT>,
		total_amount: f64,
		derivative: std::cell::Cell<StateVector::<CONSTANT>>
	}
	impl<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> FnOnce<(&[f64],)> for Derivative<Rate, CONSTANT> {
		type Output = Option<Box<[f64]>>;
		extern "rust-call" fn call_once(mut self, args: (&[f64],)) -> Self::Output { self.call_mut(args) }
	}
	impl<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> FnMut<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call_mut(&mut self, args: (&[f64],)) -> Self::Output { self.call(args) }
	}
	impl<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> Fn<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call(&self, (u,): (&[f64],)) -> Self::Output {
			let Self{rate, constant, total_amount, derivative} = self;
			get_mut(derivative, |mut derivative| {
				let pressure_R = constant.0 as f64;
				rate(*constant, &StateVector(explicit(*total_amount, pressure_R, u)), &mut derivative);
				Some(implicit(&derivative.0))
			})
		}
	}
	let ref derivative = Derivative{
		rate,
		constant,
		total_amount,
		derivative
	};

	let start = std::time::Instant::now();
	let (mut time, mut state) = (0., &*state);
	while time < 1. {
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(derivative, next_time, state); }
	}
	use itertools::Itertools;
	println!("{:.1}s {:4.1}ms {:4.0}K {} min:{:.0e}", (std::time::Instant::now()-start).as_secs_f32(), time*1e3, state[0], state[1..].iter().map(|&v| if num::abs(v)<1e-3 { if v<0. { "-0" } else { "0" }.to_owned() } else { format!("{:.1}", v) }).format(" "), state[2..].iter().copied().reduce(f64::min).unwrap());
}
