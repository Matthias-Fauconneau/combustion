#![allow(mixed_script_confusables, non_snake_case, incomplete_features, confusable_idents)]
#![feature(type_ascription, array_map, non_ascii_idents, const_generics, const_evaluatable_checked, destructuring_assignment, test)]

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
include!(concat!(env!("OUT_DIR"), "/cantera.rs"));

use {fehler::throws, error::Error};
pub use combustion::*;
//pub mod cantera;
mod reaction;
//mod transport;

#[throws] fn main() {
	let model = std::fs::read("CH4+O2.ron")?;
	let model = model::Model::new(&model)?;
	let ref simulation = Simulation::new(&model)?;
	let model = Model::new(model);
	reaction::check(model, &simulation)?;
	//transport::check(&simulation)?;
}
