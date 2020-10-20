#![feature(box_syntax)]

#[fehler::throws(anyhow::Error)] fn main() {
	use {num::real, combustion::*};
	ui::app::run(plot::Plot{
		keys: box [&["T"] as &[_], &[]],
		values: vec!(State{time: real(0.), temperature: real(0.), mass_fractions: box [real(1.)]}.into(), State{time: real(1.), temperature: real(1.), mass_fractions: box [real(0.)]}.into()),
		source: std::iter::empty()
	})?
}
