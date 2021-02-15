#![allow(mixed_script_confusables, non_snake_case, incomplete_features, confusable_idents)]
#![feature(type_ascription, array_map, non_ascii_idents, const_generics, const_evaluatable_checked, destructuring_assignment, never_type)]
use {fehler::throws, error::{Error, Result}};
pub use combustion::*;
pub mod cantera;
mod reaction;
mod transport;

//#[throws]
fn main() -> Result<!> {
	let system = default()?;
	let simulation = new(&system)?;
	reaction::check(&simulation)//?;
	//transport::check(&simulation)?;
}
