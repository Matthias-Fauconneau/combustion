#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked, destructuring_assignment)]

#[fehler::throws(Box<dyn std::error::Error>)] fn main() {
	{use trace::*; rstack_self()?; unmask_SSE_exceptions(); trace_signal_floating_point_exception();}
	let system = std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let Simulation{species_names, system, pressure_R, time_step, state, ..} = Simulation::<53>::new(&system)?;
	fn from<const S: usize>(State{temperature, amounts}: State<S>) -> Box<[Box<[f64]>]> { Box::new( [Box::new([temperature]) as Box<[_]>, Box::new(amounts)] ) as Box<[_]> }
	let mut time = 0.;
	let mut app = ui::app::App::new(ui::plot::Plot::new(Box::new( [&["T"] as &[_], &species_names] ), vec![(time, from(state))]))?;
	let total_amount = state.amounts.iter().sum();
	let mut state: [f64; 1+53-1] = state.into();
	let mut cvode = cvode::CVODE::new(move |u| system.rate/*and_jacobian*/(pressure_R, u).map(|(rate, /*jacobian*/)| rate), &state);
	app.idle = Box::new( move |plot| {
		(time, state) = cvode.step(time+time_step, &state);
		plot.values.push( (time*1e9/*ns*/, from(State::<53>::new(total_amount, &state))) );
		true
	} );
	app.run()?
}
