#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked)]

#[fehler::throws(Box<dyn std::error::Error>)] fn main() {
	{use trace::*; rstack_self()?; unmask_SSE_exceptions(); trace_signal_floating_point_exception();}
	let system = std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let Simulation{species_names, system, pressure_R, time_step, mut state, ..} = Simulation::<35>::new(&system)?;
	let ref u: [f64; 1+35-1] = state.into();
	let mut cvode = cvode::CVODE::new(move |u| system.rate/*and_jacobian*/(pressure_R, u).map(|(rate, /*jacobian*/)| rate), u);
	fn from<const S: usize>(State{temperature, amounts}: State<S>) -> Box<[Box<[f64]>]> { Box::new( [Box::new([temperature]) as Box<[_]>, Box::new(amounts)] ) as Box<[_]> }
	let mut app = ui::app::App::new(ui::plot::Plot::new(Box::new( [&["T"] as &[_], &species_names] ), vec![(0., from(state))]))?;
	app.idle = Box::new( move |plot| {
		let (time, u) = cvode.step(time_step, &state.into());
		state = State::<35>::new(state.amounts.iter().sum(), &u);
		plot.values.push( (time*1e9/*ns*/, from(state)) );
		true
	} );
	app.run()?
}
