#![allow(incomplete_features, non_snake_case)]
#![feature(const_generics, const_evaluatable_checked, destructuring_assignment, array_map, unboxed_closures, fn_traits, type_ascription, trait_alias, box_patterns)]
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

#[fehler::throws(anyhow::Error)] fn main() {
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
				rate(*constant, &StateVector(explicit(*total_amount, pressure_R, u)), &mut derivative);
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
		vec![vec![*temperature, /* *pressure_R*K*NA,*/ *volume].into_boxed_slice(), iter::map(&**amounts, |v| v/amounts.iter().sum::<f64>())].into_boxed_slice()
	}
	let state : StateVector<{Pressure}> = state.into();
	let state = implicit(&state);
	let mut cvode = cvode::CVODE::new(/*relative_tolerance:*/ 1e-4, /*absolute_tolerance:*/ 1e-7, &state);
	let plot_min_time = 0.4;
	let (mut time, mut state) = (0., &*state);
	while time < plot_min_time {
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(derivative, next_time, state); }
	}
	let mut min = 0f64;
	let values = iter::eval(1000, |_| {
		let value = ((time-plot_min_time)*1e3, from(&State::new(total_amount, constant, &StateVector::<{Pressure}>(explicit(total_amount, constant.0 as f64, &state)))));
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(derivative, next_time, state); }
		min = min.min(state[2..].iter().copied().reduce(f64::min).unwrap());
		value
	});
	println!("{:.0e}", min);
	let key = iter::map((0..len).map(|k| values.iter().map(move |(_, sets)| sets[1][k])), |Y| Y.reduce(f64::max).unwrap());
	let mut s = iter::eval(len, |i| i);
	let (_, _, select) = s.select_nth_unstable_by(len-6, |&a,&b| key[a].partial_cmp(&key[b]).unwrap());
	let species_names = iter::map(&*select, |&k| species_names[k]);
	let values = iter::map(Vec::from(values).into_iter(), |(t, sets)| {
		use std::convert::TryInto;
		let sets:Box<[_;2]> = sets.try_into().unwrap();
		let box [sets_0, sets_1] = sets;
		(t, vec![sets_0, iter::map(&*select, |&k| sets_1[k])].into_boxed_slice())
	});
	let mut plot = ui::plot::Plot::new(vec![&["T",/*"P",*/"V"] as &[_], &species_names].into_boxed_slice(), values.to_vec());
	let mut target = image::Image::zero((3840, 2160).into());
	plot.paint(&mut target.as_mut())?;
	pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
	png::save_buffer("/var/tmp/image.png", as_bytes(&target.data), target.size.x, target.size.y, png::ColorType::Rgba8).unwrap();
	//ui::app::run(plot)?
}
