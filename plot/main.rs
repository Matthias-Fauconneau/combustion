#![allow(incomplete_features, non_snake_case)]
#![feature(const_generics, const_evaluatable_checked, destructuring_assignment, array_map, unboxed_closures, fn_traits, type_ascription, trait_alias)]
fn get_mut<T: Default, U>(cell: &std::cell::Cell<T>, f: impl FnOnce(&mut T) -> U) -> U {
	let mut value = cell.take();
	let result = f(&mut value);
	cell.set(value);
	result
}

struct Row<W, const N: usize>([W; N]);
use {xy::xy, ui::widget::{size, Target, Widget, Event, EventContext}};
impl<W:Widget, const N: usize> Widget for Row<W,N> {
	fn size(&mut self, size: size) -> size { self.0[0].size(size) }
	#[fehler::throws(anyhow::Error)] fn paint(&mut self, target: &mut Target) {
			for (column, widget) in self.0.iter_mut().enumerate() {
				widget.paint(&mut target.slice_mut(xy{x: (column as u32)*target.size.x/(N as u32), y: 0}, target.size/xy{x: N as u32, y: 1}))?
			}
	}
	#[fehler::throws(anyhow::Error)] fn event(&mut self, size: size, event_context: &EventContext, event: &Event) -> bool {
		for (_column, widget) in self.0.iter_mut().enumerate() { if widget.event(size, event_context, event)? { return true; } } // todo: translate position
		false
	}
}

fn implicit(u: &[f64]) -> Box<[f64]> { [u[0]].iter().chain(&u[2..]).copied().collect() } // one of P or V imply the other using ideal gas law
fn explicit(total_amount: f64, pressure_R: f64, u: &[f64]) -> Box<[f64]> { // Reconstructs P|V using ideal gas law
	let temperature = u[0];
	let volume = total_amount / pressure_R * temperature;
	[temperature, volume].iter().chain(&u[1..]).copied().collect()
}

fn main() -> anyhow::Result<()> {
	let model = &std::fs::read("CH4+O2.ron").unwrap();
	use combustion::{*, reaction::{*, Property::*}};
	let model = model::Model::new(&model).unwrap();
	let (ref species_names, ref species) = Species::new(&model.species);
	let len = species.len();
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(species_names, r));
	let rate = rate_function(species, &*reactions);
	let ref state = combustion::initial_state(&model);
	let total_amount = state.amounts.iter().sum();

	struct Derivative<Rate: reaction::Rate<CONSTANT>, const CONSTANT: Property> {
		rate: Rate,
		constant: Constant<CONSTANT>,
		total_amount: f64,
		derivative: std::cell::Cell<StateVector::<CONSTANT>>
	}
	impl<Rate: reaction::Rate<CONSTANT>, const CONSTANT: Property> FnOnce<(&[f64],)> for Derivative<Rate, CONSTANT> {
		type Output = Option<Box<[f64]>>;
		extern "rust-call" fn call_once(mut self, args: (&[f64],)) -> Self::Output { self.call_mut(args) }
	}
	impl<Rate: reaction::Rate<CONSTANT>, const CONSTANT: Property> FnMut<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call_mut(&mut self, args: (&[f64],)) -> Self::Output { self.call(args) }
	}
	impl<Rate: reaction::Rate<CONSTANT>, const CONSTANT: Property> Fn<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call(&self, (u,): (&[f64],)) -> Self::Output {
			let Self{rate, constant, total_amount, derivative} = self;
			get_mut(derivative, |mut derivative| {
				let pressure_R = constant.0 as f64;
				let mut debug = vec![f64::NAN; /*model.reactions.len()*/325*2].into_boxed_slice();
				rate(*constant, &StateVector(explicit(*total_amount, pressure_R, u)), &mut derivative, &mut debug);
				Some(implicit(&derivative.0))
			})
		}
	}
	let constant = state.constant();
	let derivative = /*Derivative*/StateVector(std::iter::repeat(0.).take(2+len-1).collect());
	let derivative = std::cell::Cell::new(derivative);
	let ref derivative = Derivative{
		rate,
		constant,
		total_amount,
		derivative
	};

	fn from(State{temperature, pressure_R: _, volume, amounts, ..}: &State) -> Box<[Box<[f64]>]> {
		Box::new( [Box::new([*temperature, /* *pressure_R*K*NA,*/ *volume]) as Box<[_]>, iter::map(&**amounts, |v| v/amounts.iter().sum::<f64>())] ) as Box<[_]>
	}
	let state : StateVector<{Pressure}> = state.into();
	let state = implicit(&state);
	let mut cvode = cvode::CVODE::new(/*relative_tolerance:*/ 1e-4, /*absolute_tolerance:*/ 1e-7, &state);
	let (mut time, mut state) = (0., &*state);
	let mut plot = ui::plot::Plot::new(Box::new( [&["T",/*"P",*/"V"] as &[_], species_names] ), vec![]);
	let mut min = 0f64;
	while time < 0.5 {
		let plot_min_time = 0.4;
		if time > plot_min_time {
			plot.values.push(((time-plot_min_time)*1e3, from(&State::new(total_amount, constant, &StateVector::<{Pressure}>(explicit(total_amount, constant.0 as f64, &state))))));
		}
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(derivative, next_time, state); }
		min = min.min(state[2..].iter().copied().reduce(f64::min).unwrap());
	}
	println!("{:.0e}", min);
	ui::app::run(plot)
}
