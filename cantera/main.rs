#![allow(mixed_script_confusables, incomplete_features, non_snake_case)]#![feature(unboxed_closures, destructuring_assignment, fn_traits, const_generics)]
#[cfg(feature="reaction")] mod reaction;
#[cfg(feature="transport")] mod transport;
fn main() {
	#[cfg(feature="trace")] { trace::rstack_self().unwrap(); trace::signal_interrupt(); }
	#[cfg(feature="reaction")] {
		let model = std::fs::read("CH4.ron").unwrap();
		let ref model = combustion::model::Model::new(&model).unwrap();
		let ref state = combustion::initial_state(model);
		reaction::check(model, state);
	}
	#[cfg(feature="transport")] transport::check("/usr/share/cantera/data/LiDryer.yaml");
}
