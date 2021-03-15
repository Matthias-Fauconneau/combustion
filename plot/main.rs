#![allow(incomplete_features, non_snake_case)]
#![feature(const_generics, const_evaluatable_checked, destructuring_assignment, array_map, unboxed_closures, fn_traits, type_ascription)]
use {fehler::throws, error::Error, combustion::{*, Property::*}};

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

#[throws] fn main() {
	let model = &std::fs::read("CH4+O2.ron")?;
	let model = model::Model::new(&model)?;
	let Simulation{ref species_names, time_step, state} = Simulation::new(&model)?;
	let model = Model::new(model);
	let len = model.len();
	let (_traps, (_function, _size), rate) = model.rate();
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
	//let mut state: StateVector<{Pressure}> = (&state).into();
	let mut state = ideal(&promote(&((&state).into():StateVector<{Pressure}>)));
	pub fn map<T, U, F: Fn(&T)->U>(v: &[T], f: F) -> Box<[U]> { v.iter().map(f).collect() }
	pub fn promote(v: &[f32]) -> Box<[f64]> { map(v, |&v| v as f64) }
	pub fn demote(v: &[f64]) -> Box<[f32]> { map(v, |&v| v as f32) }
	fn ideal(u: &[f64]) -> Box<[f64]> { [u[0]].iter().chain(&u[2..]).copied().collect() } // Skip extra state variable (P|V): for ideal gas (PV=nRT) one constant directly imply the other through T
	fn real(total_amount: f64, pressure: f64, u: &[f64]) -> Box<[f64]> { // Reconstructs extra state variable (P|V) using ideal gas law
		let temperature = u[0];
		let volume = K * total_amount / pressure * temperature;
		[temperature, volume].iter().chain(&u[1..]).copied().collect()
	}
	let mut cvode = cvode::CVODE::new(&state);
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
			let mut derivative = self.derivative.take();
			let pressure = self.constant.0 as f64;
			(self.rate)(self.constant, &StateVector(map(&real(self.total_amount, pressure, u), |&v| (v as f32).max(0.))), &mut derivative);
			let result = Some(ideal(&promote(&derivative.0)));
			self.derivative.set(derivative);
			//eprintln!("{:e}", result.as_ref().unwrap()[0]);
			result
		}
	}
	let derivative = /*Derivative*/StateVector::<{Pressure}>(std::iter::repeat(0.).take(2+len-1).collect());
	let derivative = Derivative{
		rate,
		constant,
		total_amount,
		derivative: std::cell::Cell::new(derivative)
	};
	app.run(|plots| { //: &mut Row<[ui::plot::Plot; 2]>
		for _ in 0..100 { (time, state) = {
			let (time, state) = cvode.step(&derivative, time+time_step, &state);
			(time, state.to_vec().into_boxed_slice())
		}}
		plots.0[0].values.push((time*1e3/*ms*/, from(&State::new(total_amount, constant, &StateVector(demote(&real(total_amount, constant.0 as f64, &state)))))));
		//plots.0[1].values.push((time*1e9/*ns*/, from(&reaction(simulation, time)?)));
		Ok(true)
		//Err(error::anyhow!(""))
	})?
}
