#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	use combustion::*;
	let Simulation{species: _, system, mut state} = Simulation::new(&system)?;
	/*while state.time < num::real(4e-6) {
		for _ in 0..10000 { state.step(&system); }
	}*/
	state.step(&system);
}
