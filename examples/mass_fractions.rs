#![feature(non_ascii_idents)] #![allow(mixed_script_confusables)] #[fehler::throws(anyhow::Error)] fn main() {
	use {itertools::Itertools, combustion::*};
	let system = std::fs::read("H2+O2.ron")?;
	let ron::System{species, reactions, phases, time_step} = ::ron::de::from_bytes(&system)?;
	let specie_names = collect(species.keys().copied());
	let molar_masses = collect(species.values().map(|Specie{composition,..}| composition.iter().map(|(element, &count)| (count as f64)*standard_atomic_weights[element]).sum()));
	let masses = molar_masses.iter().zip(&amounts).map(|(w,n)| w*n).collect::<Box<_>>();
	let mass = masses.iter().sum::<f64>();
	let mass_fractions = masses.iter().map(|m| m/mass);
	println!("{:?}", species.iter().zip(mass_fractions).format(" "));
}
