#![feature(type_ascription, array_map, non_ascii_idents)]#![allow(mixed_script_confusables,non_snake_case)]
extern "C" {
	fn cantera(relative_tolerance: f64, absolute_tolerance: f64, temperature: &mut f64, pressure: &mut f64, mole_proportions: *const std::os::raw::c_char, time_step: f64,
										species_len: &mut usize, species: &mut *const *const std::os::raw::c_char, concentrations: &mut *const f64);
}

#[fehler::throws(anyhow::Error)] fn main() {
	use combustion::*;
	let system = std::fs::read("H2+O2.ron")?;
	pub const S : usize = 9; // Total number of species
	pub const N : usize = 2/*T, V*/+S-1; // Skips most abundant specie (last index) (will be deduced from conservation)
	let Simulation{species: species_name, system, state: State{temperature, ref amounts}, volume, pressure, time_step, ..} = Simulation::<S,{S-1},N>::new(&system)?;
	let mut state : [_; N] = {
			use {std::convert::TryInto, iter::{into::IntoChain, array_from_iter as from_iter}};
			from_iter([temperature,volume].chain(amounts[..S-1].try_into().unwrap():[_;S-1]))
	};
	use itertools::Itertools;
	println!("{}", species_name.iter().map(|k| format!("{:^15}",k)).format(" "));
	for _ in 0..2 {
		let (other_T, ref other) = {
			let (T, _, ref concentrations) : (_,_,[_;S]) = System::<S,{S-1},N>::state(pressure, &state);
			let mole_proportions = format!("{}", species_name.iter().zip(concentrations).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
			let mole_proportions = std::ffi::CString::new(mole_proportions).unwrap();
			use std::ptr::null;
			let (mut T, mut pressure, mut species_len, mut species, mut concentrations) = (T, pressure, 0, null(), null());
			unsafe {
				cantera(/*relative_tolerance:*/ 1e-6, /*absolute_tolerance:*/ 1e-10, &mut T, &mut pressure, mole_proportions.as_ptr(), time_step, &mut species_len, &mut species, &mut concentrations);
				let species = iter::box_collect(std::slice::from_raw_parts(species, species_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
				let order = |o:Box<[_]>| iter::vec::eval(species_name, |s| o[species.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}",s,species))]);
				let concentrations = iter::box_collect(std::slice::from_raw_parts(concentrations, species_len).iter().map(|c| c*1000.)); // kmol/m^3 => mol/m^3
				(T, order(concentrations))
			}
		};
		state = system.step(/*relative_tolerance:*/ 1e-6, /*absolute_tolerance:*/ 1e-10, time_step, pressure, state);
		let (T, _, ref concentrations) : (_,_,[_;S]) = System::<S,{S-1},N>::state(pressure, &state);
		println!("{:.5} {:.5}", T, other_T);
		println!("{}", &[concentrations, other].iter().format_with("\n",|concentrations, f| f(&format_args!("{}",concentrations.iter().map(|c| format!("{:15.6e}",c)).format(" ")))));
		println!("{}", concentrations.iter().zip(other).filter(|&(&c,_)| c > 0.).map(|(c,o)| format!("{:15.8e}",c-o)).format(" "));
		use num::abs;
		println!("{}", concentrations.iter().zip(other).filter(|&(&c,_)| c > 0.).map(|(&c,&o)| format!("{:15.9e}",abs(c-o)/c.min(o))).format(" "));
	}
}
