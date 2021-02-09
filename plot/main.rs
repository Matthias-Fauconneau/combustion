#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked, destructuring_assignment, array_map, bindings_after_at)]
use {fehler::throws, error::Error, combustion::*};

mod cantera {
extern "C" {
pub fn reaction(pressure: &mut f64, temperature: &mut f64, mole_proportions: *const std::os::raw::c_char,
													species_len: &mut usize, species: &mut *const *const std::os::raw::c_char,
													rate_time: f64,
														net_productions_rates: &mut *const f64,
														reactions_len: &mut usize,
															equations: &mut *const *const std::os::raw::c_char,
															equilibrium_constants: &mut *const f64,
															forward: &mut *const f64,
															reverse: &mut *const f64,
													state_time: f64,
													concentrations: &mut *const f64);
}
}

#[throws] pub fn reaction<const S: usize>(Simulation{species_names, pressure_R, state, ..}: &Simulation<S>, time: f64) -> State<S>
where [(); S-1]:, [(); 1+S-1]: {
	use itertools::Itertools;
	let mole_proportions = format!("{}", species_names.iter().zip(&state.amounts).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
	let mole_proportions = std::ffi::CString::new(mole_proportions)?;
	use std::ptr::null;
	let mut pressure = pressure_R * (combustion::kB*combustion::NA);
	let mut temperature = state.temperature;
	let (mut species_len, mut specie_names, mut net_productions_rates, mut concentrations,
				mut reactions_len, mut equations, mut equilibrium_constants, [mut forward, mut reverse]) = (0, null(), null(), null(), 0, null(), null(), [null(); 2]);
	unsafe {
		cantera::reaction(&mut pressure, &mut temperature, mole_proportions.as_ptr(), &mut species_len, &mut specie_names,
																time, &mut net_productions_rates, &mut reactions_len, &mut equations, &mut equilibrium_constants, &mut forward, &mut reverse,
																time, &mut concentrations);
		let specie_names = iter::box_collect(std::slice::from_raw_parts(specie_names, species_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
		let order = |o:&[_]| iter::vec::eval(species_names, |s| o[specie_names.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}", s, species_names))]);
		let concentrations = order(std::slice::from_raw_parts(concentrations, species_len)).map(|c| c*1000.); // kmol/m^3 => mol/m^3
		State{temperature, amounts: iter::vec::eval(concentrations, |c| c * System::<S>::volume)}
	}
}

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
	{use trace::*; rstack_self()?; unmask_SSE_exceptions(); trace_signal_floating_point_exception();}
	let system = std::fs::read("CH4+O2.ron")?;
	let ref _simulation@Simulation{ref species_names, ref system, ref state, ref pressure_R, ref time_step, ..} = Simulation::<53>::new(&system)?;
	fn from<const S: usize>(State{temperature, amounts}: &State<S>) -> Box<[Box<[f64]>]> {
		Box::new( [Box::new([*temperature]) as Box<[_]>, Box::new(*amounts)] ) as Box<[_]>
	}
	let mut time = 0.;
	let total_amount = state.amounts.iter().sum();
	let plot = || ui::plot::Plot::new(Box::new( [&["T"] as &[_], species_names] ), vec![(time, from(state))]);
	let app = ui::app::App::new(Row([plot()/*, plot()*/]))?;
	let mut state: [f64; 1+53-1] = (*state).into();
	let mut cvode = cvode::CVODE::new(&state);
	app.run(|plots| { //: &mut Row<[ui::plot::Plot; 2]>
		for _ in 0..1000 {
			let temperature = state[0];
			let C0 = pressure_R / temperature;
			(time, state) = cvode.step(move |u| system.rate/*and_jacobian*/(*pressure_R, u).map(|(rate, /*jacobian*/)| rate), time+time_step, &state);
			// Constant pressure, Constant volume => flow on temperature change
			let temperature = state[0];
			let C = pressure_R / temperature;
			for n in &mut state[1..] { *n *= C/C0; }
		}
		plots.0[0].values.push((time*1e9/*ns*/, from(&State::<53>::new(total_amount, &state))));
		//plots.0[1].values.push((time*1e9/*ns*/, from(&reaction(simulation, time)?)));
		Ok(true)
	})?
}
