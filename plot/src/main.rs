#![feature(box_syntax)]

#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	use combustion::*;
	let Simulation{species, system, mut state} = Simulation::new(&system)?;
	let mut app = ui::app::App::new(plot::Plot::new(box [&["T"] as &[_], &species], vec!(state.clone().into())))?;
	app.idle = box move |plot| {
		for _ in 0..10000 { state.step(&system); }
		plot.values.push(state.clone().into());
		state.time < num::real(4e-6)
	};
	app.run()?
}
