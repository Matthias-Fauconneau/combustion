#![allow(incomplete_features, non_snake_case)]
#![feature(const_generics, const_evaluatable_checked, destructuring_assignment, array_map, unboxed_closures, fn_traits, type_ascription)]
use {fehler::throws, error::Error, combustion::{*, reaction::{*, Property::*}}};

struct Row<W, const N: usize>([W; N]);
use {xy::xy, ui::widget::{size, Target, Widget, Event, EventContext}};
impl<W:Widget, const N: usize> Widget for Row<W,N> {
	fn size(&mut self, size: size) -> size { self.0[0].size(size) }
	#[throws] fn paint(&mut self, target: &mut Target) {
			for (column, widget) in self.0.iter_mut().enumerate() {
				widget.paint(&mut target.slice_mut(xy{x: (column as u32)*target.size.x/(N as u32), y: 0}, target.size/xy{x: N as u32, y: 1}))?
			}
	}
	#[throws] fn event(&mut self, size: size, event_context: &EventContext, event: &Event) -> bool {
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

#[throws] fn main() {
	let model = &std::fs::read("CH4+O2.ron")?;
	let model = model::Model::new(&model)?;
	let (ref species_names, ref species) = Species::new(&model.species);
	let len = species.len();
	let (_, rate) = rate(species, iter::map(&model.reactions, |r| Reaction::new(species_names, r)));
	let ref state = combustion::initial_state(model);
	fn from(State{/*temperature, pressure, volume,*/ amounts, ..}: &State) -> Box<[Box<[f64]>]> {
		//eprintln!("{} {} {}", temperature, *pressure*NA, *volume);
		Box::new( [/*Box::new([*temperature, *pressure*NA, *volume]) as Box<[_]>,*/ amounts.clone()] ) as Box<[_]>
	}
	let mut time = 0.;
	let total_amount = state.amounts.iter().sum();
	//let plot = || ui::plot::Plot::new(Box::new( [&["T","P","V"] as &[_], species_names] ), vec![(time, from(&state))]);
	let plot = || ui::plot::Plot::new(Box::new( [species_names] ), vec![(time, from(&state))]);
	let app = ui::app::App::new(Row([plot()/*, plot()*/]))?;
	let constant = state.constant();
	let mut state = implicit(&(state.into():StateVector<{Volume}>));
	let mut cvode = cvode::CVODE::new(&state);
	let derivative = /*Derivative*/StateVector(std::iter::repeat(0.).take(2+len-1).collect());
	let derivative = std::cell::Cell::new(derivative);
	let derivative = move |u| {
		let mut derivative = derivative.take();
		let pressure = constant.0 as f64;
		rate(constant, &StateVector(iter::map(&explicit(total_amount, pressure, u), |&v| (v as f32).max(0.))), &mut derivative);
		let result = Some(implicit(&derivative.0));
		derivative.set(derivative);
		//eprintln!("{:e}", result.as_ref().unwrap()[0]);
		result
	};
	app.run(move |plots| { //: &mut Row<[ui::plot::Plot; 2]>
		for _ in 0..100 { (time, state) = {
			let (time, state) = cvode.step(&derivative, time+model.time_step, &state);
			(time, state.to_vec().into_boxed_slice())
		}}
		plots.0[0].values.push((time*1e3/*ms*/, from(&State::new(total_amount, constant, &StateVector(explicit(total_amount, constant.0 as f64, &state))))));
		//plots.0[1].values.push((time*1e9/*ns*/, from(&reaction(simulation, time)?)));
		Ok(true)
		//Err(error::anyhow!(""))
	})?
}
