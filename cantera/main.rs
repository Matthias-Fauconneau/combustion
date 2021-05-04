#![allow(mixed_script_confusables, non_snake_case, incomplete_features, confusable_idents, uncommon_codepoints)]
#![feature(type_ascription, array_map, const_generics, const_evaluatable_checked, destructuring_assignment, test, unboxed_closures, fn_traits)]

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
//include!(concat!(env!("OUT_DIR"), "/cantera.rs"));

use {fehler::throws, error::Error};
#[cfg(feature="reaction")] mod reaction;
#[cfg(feature="transport")] mod transport;

#[throws] fn main() {
	trace::rstack_self()?; trace::signal_interrupt();
	let model = std::fs::read("CH4+O2.ron")?;
	let ref model = combustion::model::Model::new(&model)?;
	let ref state = combustion::initial_state(model);
	#[cfg(feature="reaction")] reaction::check(model, state)?;
	#[cfg(feature="transport")] transport::check(model, state);
}
