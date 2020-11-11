#![feature(non_ascii_idents)] #![allow(mixed_script_confusables)] #[fehler::throws(anyhow::Error)] fn main() {
	use {itertools::Itertools, combustion::*};
	let system = std::fs::read("H2+O2.ron")?;
	pub const S : usize = 9; // Total number of species
	pub const N : usize = 2/*T, V*/+S-1; // Skips most abundant specie (last index) (will be deduced from conservation)
	let Simulation{species, system: System{molar_masses, ..}, state: State{amounts, ..}, ..} = Simulation::<S,{S-1},N>::new(&system)?;
	let masses = molar_masses.iter().zip(&amounts).map(|(w,n)| w*n).collect::<Box<_>>();
	let mass = masses.iter().sum::<f64>();
	let mass_fractions = masses.iter().map(|m| m/mass);
	println!("{:?}", species.iter().zip(mass_fractions).format(" "));
}
