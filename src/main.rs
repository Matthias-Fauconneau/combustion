#![feature(box_syntax)]

#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	use combustion::*;
	let Simulation{species, system, mut state} = Simulation::new(&system)?;
	let mut app = ui::app::App::new(plot::Plot::new(box [&["T"] as &[_], &species], vec!(state.clone().into())))?;
	app.idle = box move |plot| {
		state.step(&system);
		plot.values.push(state.clone().into());
		true
	};
	app.run()?
}
